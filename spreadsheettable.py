# -*- coding: utf-8 -*-

__author__ = u'Tomasz Świderski <contact@tomaszswiderski.com>'
__copyright__ = u'Copyright (c) 2010 Tomasz Świderski'

from reportlab.platypus.tables import *
from reportlab.platypus.tables import (_rowLen, _calc_pc, _hLine, _multiLine,
    _convert2int, _endswith, _isLineCommand, _setCellStyle)

from formula import Formula

def spanFixDim(V0,V,spanCons,FUZZ=rl_config._FUZZ):
    #assign required space to variable rows equally to existing calculated values
    M = {}
    for (x0,x1),v in spanCons.iteritems():
        t = sum([V[x]+M.get(x,0) for x in xrange(x0,x1+1)])
        if t>=v-FUZZ: continue      #already good enough
        X = [x for x in xrange(x0,x1+1) if V0[x] is None]   #variable candidates
        if not X: continue          #something wrong here mate
        v -= t
        v /= float(len(X))
        for x in X:
            M[x] = M.get(x,0) + v
    for x,v in M.iteritems():
        V[x] += v


class SpreadsheetTable(Flowable):
    def __init__(self, data, colWidths=None, rowHeights=None, style=None,
        repeatRows=0, repeatCols=0, splitByRow=1, emptyTableAction=None,
        ident=None, hAlign='CENTER', vAlign='MIDDLE', normalizedData=0,
        cellStyles=None, activeRows=None, repeatRowsB=0):

        self.ident = ident
        self.hAlign = hAlign
        self.vAlign = vAlign
        if not isinstance(data,(tuple,list)):
            raise ValueError("%s invalid data type" % self.identity())
        self._nrows = nrows = len(data)
        self.repeatRows = repeatRows
        self._repeatRowsB = repeatRowsB
        if activeRows is None:
            activeRows = (None, nrows - repeatRowsB)
        self._activeRows = activeRows
        self._cellvalues = []
        _seqCW = isinstance(colWidths,(tuple,list))
        _seqRH = isinstance(rowHeights,(tuple,list))
        if nrows: self._ncols = ncols = max(map(_rowLen,data))
        elif colWidths and _seqCW: ncols = len(colWidths)
        else: ncols = 0
        if not emptyTableAction: emptyTableAction = rl_config.emptyTableAction
        if not (nrows and ncols):
            if emptyTableAction=='error':
                raise ValueError("%s must have at least a row and column" % self.identity())
            elif emptyTableAction=='indicate':
                self.__class__ = Preformatted
                global _emptyTableStyle
                if '_emptyTableStyle' not in globals().keys():
                    _emptyTableStyle = ParagraphStyle('_emptyTableStyle')
                    _emptyTableStyle.textColor = colors.red
                    _emptyTableStyle.backColor = colors.yellow
                Preformatted.__init__(self,'%s(%d,%d)' % (self.__class__.__name__,nrows,ncols), _emptyTableStyle)
            elif emptyTableAction=='ignore':
                self.__class__ = Spacer
                Spacer.__init__(self,0,0)
            else:
                raise ValueError('%s bad emptyTableAction: "%s"' % (self.identity(),emptyTableAction))
            return

        # we need a cleanup pass to ensure data is strings - non-unicode and non-null
        if normalizedData:
            self._cellvalues = data
        else:
            self._cellvalues = data = self.normalizeData(data)
        if not _seqCW: colWidths = ncols*[colWidths]
        elif len(colWidths)!=ncols:
            if rl_config.allowShortTableRows and isinstance(colWidths,list):
                n = len(colWidths)
                if n<ncols:
                    colWidths[n:] = (ncols-n)*[colWidths[-1]]
                else:
                    colWidths = colWidths[:ncols]
            else:
                raise ValueError("%s data error - %d columns in data but %d in column widths" % (self.identity(),ncols, len(colWidths)))
        if not _seqRH: rowHeights = nrows*[rowHeights]
        elif len(rowHeights) != nrows:
            raise ValueError("%s data error - %d rows in data but %d in row heights" % (self.identity(),nrows, len(rowHeights)))
        for i,d in enumerate(data):
            n = len(d)
            if n!=ncols:
                if rl_config.allowShortTableRows and isinstance(d,list):
                    d[n:] = (ncols-n)*['']
                else:
                    raise ValueError("%s expected %d not %d columns in row %d!" % (self.identity(),ncols,n,i))
        self._argH = self._rowHeights = rowHeights
        self._colWidths = self._argW = colWidths
        if cellStyles is None:
            cellrows = []
            for i in xrange(nrows):
                cellcols = []
                for j in xrange(ncols):
                    cellcols.append(CellStyle(`(i,j)`))
                cellrows.append(cellcols)
            self._cellStyles = cellrows
        else:
            self._cellStyles = cellStyles

        self._bkgrndcmds = []
        self._linecmds = []
        self._spanCmds = []
        self._nosplitCmds = []
        self.repeatCols = repeatCols
        self.splitByRow = splitByRow

        if style:
            self.setStyle(style)

    def __repr__(self):
        "incomplete, but better than nothing"
        r = getattr(self,'_rowHeights','[unknown]')
        c = getattr(self,'_colWidths','[unknown]')
        cv = getattr(self,'_cellvalues','[unknown]')
        import pprint
        cv = pprint.pformat(cv)
        cv = cv.replace("\n", "\n  ")
        return "%s(\n rowHeights=%s,\n colWidths=%s,\n%s\n) # end table" % (self.__class__.__name__,r,c,cv)

    def normalizeData(self, data):
        """Takes a block of input data (list of lists etc.) and
        - coerces unicode strings to non-unicode UTF8
        - coerces nulls to ''
        """
        def normCell(stuff):
            if stuff is None:
                return ''
            elif isinstance(stuff,unicode):
                return stuff.encode('utf8')
            else:
                return stuff
        outData = []
        for row in data:
            outRow = [normCell(cell) for cell in row]
            outData.append(outRow)
        return outData

    def identity(self, maxLen=30):
        '''Identify our selves as well as possible'''
        if self.ident: return self.ident
        vx = None
        nr = getattr(self,'_nrows','unknown')
        nc = getattr(self,'_ncols','unknown')
        cv = getattr(self,'_cellvalues',None)
        if cv and 'unknown' not in (nr,nc):
            b = 0
            for i in xrange(nr):
                for j in xrange(nc):
                    v = cv[i][j]
                    if isinstance(v,(list,tuple,Flowable)):
                        if not isinstance(v,(tuple,list)): v = (v,)
                        r = ''
                        for vij in v:
                            r = vij.identity(maxLen)
                            if r and r[-4:]!='>...':
                                break
                        if r and r[-4:]!='>...':
                            ix, jx, vx, b = i, j, r, 1
                    else:
                        v = v is None and '' or str(v)
                        ix, jx, vx = i, j, v
                        b = (vx and isinstance(v,basestring)) and 1 or 0
                        if maxLen: vx = vx[:maxLen]
                    if b: break
                if b: break
        if vx:
            vx = ' with cell(%d,%d) containing\n%s' % (ix,jx,repr(vx))
        else:
            vx = '...'

        return "<%s@0x%8.8X %s rows x %s cols>%s" % (self.__class__.__name__, id(self), nr, nc, vx)

    def _listCellGeom(self, V,w,s,W=None,H=None,aH=72000):
        if not V: return 0,0
        aW = w - s.leftPadding - s.rightPadding
        aH = aH - s.topPadding - s.bottomPadding
        t = 0
        w = 0
        canv = getattr(self,'canv',None)
        sb0 = None
        for v in V:
            vw, vh = v.wrapOn(canv, aW, aH)
            sb = v.getSpaceBefore()
            sa = v.getSpaceAfter()
            if W is not None: W.append(vw)
            if H is not None: H.append(vh)
            w = max(w,vw)
            t += vh + sa + sb
            if sb0 is None:
                sb0 = sb
        return w, t - sb0 - sa

    def _listValueWidth(self,V,aH=72000,aW=72000):
        if not V: return 0,0
        t = 0
        w = 0
        canv = getattr(self,'canv',None)
        return max([v.wrapOn(canv,aW,aH)[0] for v in V])

    def _calc_width(self,availWidth,W=None):
        if getattr(self,'_width_calculated_once',None): return
        #comments added by Andy to Robin's slightly terse variable names
        if not W: W = _calc_pc(self._argW,availWidth)   #widths array
        if None in W:  #some column widths are not given
            canv = getattr(self,'canv',None)
            saved = None
            if self._spanCmds:
                colSpanCells = self._colSpanCells
                spanRanges = self._spanRanges
            else:
                colSpanCells = ()
                spanRanges = {}
            spanCons = {}
            if W is self._argW:
                W0 = W
                W = W[:]
            else:
                W0 = W[:]
            V = self._cellvalues
            S = self._cellStyles
            while None in W:
                j = W.index(None) #find first unspecified column
                w = 0
                for i,Vi in enumerate(V):
                    v = Vi[j]
                    s = S[i][j]
                    ji = j,i
                    span = spanRanges.get(ji,None)
                    if ji in colSpanCells and not span: #if the current cell is part of a spanned region,
                        t = 0.0                         #assume a zero size.
                    else:#work out size
                        t = self._elementWidth(v, s, (j, i))
                        if t is None:
                            raise ValueError("Flowable %s in cell(%d,%d) can't have auto width\n%s" % (v.identity(30),i,j,self.identity(30)))
                        t += s.leftPadding+s.rightPadding
                        if span:
                            c0 = span[0]
                            c1 = span[2]
                            if c0!=c1:
                                x = c0,c1
                                spanCons[x] = max(spanCons.get(x,t),t)
                                t = 0
                    if t>w: w = t   #record a new maximum

                W[j] = w

            if spanCons:
                spanFixDim(W0,W,spanCons)

        self._colWidths = W
        width = 0
        self._colpositions = [0]        #index -1 is right side boundary; we skip when processing cells
        for w in W:
            width = width + w
            self._colpositions.append(width)

        self._width = width
        self._width_calculated_once = 1

    def _elementWidth(self, v, s, cellcoord):
        if isinstance(v,(list,tuple)):
            w = 0
            for e in v:
                ew = self._elementWidth(e,s)
                if ew is None: return None
                w = max(w,ew)
            return w
        elif isinstance(v,Flowable) and v._fixedWidth:
            if hasattr(v, 'width') and isinstance(v.width,(int,float)): return v.width
            if hasattr(v, 'drawWidth') and isinstance(v.drawWidth,(int,float)): return v.drawWidth
        elif isinstance(v, Formula):
            v = v.get_max_value(self._cellvalues, self.repeatRows,
                self._repeatRowsB, cellcoord)
        # Even if something is fixedWidth, the attribute to check is not
        # necessarily consistent (cf. Image.drawWidth).  Therefore, we'll
        # be extra-careful and fall through to this code if necessary.
        if hasattr(v, 'minWidth'):
            try:
                w = v.minWidth() # should be all flowables
                if isinstance(w,(float,int)): return w
            except AttributeError:
                pass
        v = (v is not None and str(v) or '').split("\n")
        fontName = s.fontname
        fontSize = s.fontsize
        return max([stringWidth(x,fontName,fontSize) for x in v])

    def _calc_height(self, availHeight, availWidth, H=None):
        H0 = self._argH
        H = self._rowHeights
        W = self._colWidths

        hmax = lim = len(H)

        if None in H:
            canv = getattr(self,'canv',None)
            saved = None
            #get a handy list of any cells which span rows. should be ignored for sizing
            if self._spanCmds:
                rowSpanCells = self._rowSpanCells
                colSpanCells = self._colSpanCells
                spanRanges = self._spanRanges
                colpositions = self._colpositions
            else:
                rowSpanCells = colSpanCells = ()
                spanRanges = {}
            if canv: saved = canv._fontname, canv._fontsize, canv._leading
            spanCons = {}
            FUZZ = rl_config._FUZZ
            while None in H:
                i = H.index(None)
                V = self._cellvalues[i] # values for row i
                S = self._cellStyles[i] # styles for row i
                h = 0
                j = 0
                for j,(v, s, w) in enumerate(zip(V, S, W)): # value, style, width (lengths must match)
                    ji = j,i
                    span = spanRanges.get(ji,None)
                    if ji in rowSpanCells and not span:
                        continue # don't count it, it's either occluded or unreliable

                    if isinstance(v,(tuple,list,Flowable)):
                        if isinstance(v,Flowable): v = (v,)
                        if w is None and not self._canGetWidth(v):
                            raise ValueError("Flowable %s in cell(%d,%d) can't have auto width in\n%s" % (v[0].identity(30),i,j,self.identity(30)))
                        if canv: canv._fontname, canv._fontsize, canv._leading = s.fontname, s.fontsize, s.leading or 1.2*s.fontsize
                        if ji in colSpanCells:
                            if not span: continue
                            w = max(colpositions[span[2]+1]-colpositions[span[0]],w)
                        dW,t = self._listCellGeom(v,w or self._listValueWidth(v),s)
                        if canv: canv._fontname, canv._fontsize, canv._leading = saved
                        dW = dW + s.leftPadding + s.rightPadding
                        if not rl_config.allowTableBoundsErrors and dW>w:
                            from reportlab.platypus.doctemplate import LayoutError
                            raise LayoutError("Flowable %s (%sx%s points) too wide for cell(%d,%d) (%sx* points) in\n%s" % (v[0].identity(30),fp_str(dW),fp_str(t),i,j, fp_str(w), self.identity(30)))
                    else:
                        v = (v is not None and str(v) or '').split("\n")
                        t = (s.leading or 1.2*s.fontSize)*len(v)
                    t += s.bottomPadding+s.topPadding
                    if span:
                        r0 = span[1]
                        r1 = span[3]
                        if r0!=r1:
                            x = r0,r1
                            spanCons[x] = max(spanCons.get(x,t),t)
                            t = 0
                    if t>h: h = t   #record a new maximum
                H[i] = h

            if spanCons:
                spanFixDim(H0,H,spanCons)

        hmax = self._activeRows[1]
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        height = self._height = sum(H[:self.repeatRows] +
            H[activeRows0:hmax] +
            H[self._nrows-self._repeatRowsB:])
        self._rowpositions = [height]    # index 0 is actually topline; we skip when processing cells
        for h in H[:self.repeatRows] + H[activeRows0:hmax] + H[self._nrows-self._repeatRowsB:]:
            height = height - h
            self._rowpositions.append(height)
        assert abs(height)<1e-8, 'Internal height error'

    def _calc(self, availWidth, availHeight):
        #if hasattr(self,'_width'): return

        #in some cases there are unsizable things in
        #cells.  If so, apply a different algorithm
        #and assign some withs in a less (thanks to Gary Poster) dumb way.
        #this CHANGES the widths array.
        if (None in self._colWidths or '*' in self._colWidths) and self._hasVariWidthElements():
            W = self._calcPreliminaryWidths(availWidth) #widths
        else:
            W = None

        # need to know which cells are part of spanned
        # ranges, so _calc_height and _calc_width can ignore them
        # in sizing
        if self._spanCmds:
            self._calcSpanRanges()

        if self._nosplitCmds:
            self._calcNoSplitRanges()

        # calculate the full table width
        self._calc_width(availWidth,W=W)

        # calculate the full table height
        self._calc_height(availHeight,availWidth)

        if self._spanCmds:
            #now work out the actual rect for each spanned cell from the underlying grid
            self._calcSpanRects()

    def _hasVariWidthElements(self, upToRow=None):
        """Check for flowables in table cells and warn up front.

        Allow a couple which we know are fixed size such as
        images and graphics."""
        if upToRow is None: upToRow = self._nrows
        for row in xrange(min(self._nrows, upToRow)):
            for col in xrange(self._ncols):
                value = self._cellvalues[row][col]
                if not self._canGetWidth(value):
                    return 1
        return 0

    def _canGetWidth(self, thing):
        "Can we work out the width quickly?"
        if isinstance(thing,(list, tuple)):
            for elem in thing:
                if not self._canGetWidth(elem):
                    return 0
            return 1
        elif isinstance(thing, Flowable):
            return thing._fixedWidth  # must loosen this up
        else: #str, number, None etc.
            #anything else gets passed to str(...)
            # so should be sizable
            return 1

    def _calcPreliminaryWidths(self, availWidth):
        """Fallback algorithm for when main one fails.

        Where exact width info not given but things like
        paragraphs might be present, do a preliminary scan
        and assign some best-guess values."""

        W = list(self._argW) # _calc_pc(self._argW,availWidth)
        verbose = 0
        totalDefined = 0.0
        percentDefined = 0
        percentTotal = 0
        numberUndefined = 0
        numberGreedyUndefined = 0
        for w in W:
            if w is None:
                numberUndefined += 1
            elif w == '*':
                numberUndefined += 1
                numberGreedyUndefined += 1
            elif _endswith(w,'%'):
                percentDefined += 1
                percentTotal += float(w[:-1])
            else:
                assert isinstance(w,(int,float))
                totalDefined = totalDefined + w
        if verbose: print 'prelim width calculation.  %d columns, %d undefined width, %0.2f units remain' % (
            self._ncols, numberUndefined, availWidth - totalDefined)

        #check columnwise in each None column to see if they are sizable.
        given = []
        sizeable = []
        unsizeable = []
        minimums = {}
        totalMinimum = 0
        elementWidth = self._elementWidth
        for colNo in xrange(self._ncols):
            w = W[colNo]
            if w is None or w=='*' or _endswith(w,'%'):
                siz = 1
                current = final = None
                for rowNo in xrange(self._nrows):
                    value = self._cellvalues[rowNo][colNo]
                    style = self._cellStyles[rowNo][colNo]
                    new = elementWidth(value,style, (colNo, rowNo))
                    new += style.leftPadding + style.rightPadding
                    final = max(current, new)
                    current = new
                    siz = siz and self._canGetWidth(value) # irrelevant now?
                if siz:
                    sizeable.append(colNo)
                else:
                    unsizeable.append(colNo)
                minimums[colNo] = final
                totalMinimum += final
            else:
                given.append(colNo)
        if len(given) == self._ncols:
            return
        if verbose: print 'predefined width:   ',given
        if verbose: print 'uncomputable width: ',unsizeable
        if verbose: print 'computable width:   ',sizeable

        # how much width is left:
        remaining = availWidth - (totalMinimum + totalDefined)
        if remaining > 0:
            # we have some room left; fill it.
            definedPercentage = (totalDefined/availWidth)*100
            percentTotal += definedPercentage
            if numberUndefined and percentTotal < 100:
                undefined = numberGreedyUndefined or numberUndefined
                defaultWeight = (100-percentTotal)/undefined
                percentTotal = 100
                defaultDesired = (defaultWeight/percentTotal)*availWidth
            else:
                defaultWeight = defaultDesired = 1
            # we now calculate how wide each column wanted to be, and then
            # proportionately shrink that down to fit the remaining available
            # space.  A column may not shrink less than its minimum width,
            # however, which makes this a bit more complicated.
            desiredWidths = []
            totalDesired = 0
            effectiveRemaining = remaining
            for colNo, minimum in minimums.items():
                w = W[colNo]
                if _endswith(w,'%'):
                    desired = (float(w[:-1])/percentTotal)*availWidth
                elif w == '*':
                    desired = defaultDesired
                else:
                    desired = not numberGreedyUndefined and defaultDesired or 1
                if desired <= minimum:
                    W[colNo] = minimum
                else:
                    desiredWidths.append(
                        (desired-minimum, minimum, desired, colNo))
                    totalDesired += desired
                    effectiveRemaining += minimum
            if desiredWidths: # else we're done
                # let's say we have two variable columns.  One wanted
                # 88 points, and one wanted 264 points.  The first has a
                # minWidth of 66, and the second of 55.  We have 71 points
                # to divide up in addition to the totalMinimum (i.e.,
                # remaining==71).  Our algorithm tries to keep the proportion
                # of these variable columns.
                #
                # To do this, we add up the minimum widths of the variable
                # columns and the remaining width.  That's 192.  We add up the
                # totalDesired width.  That's 352.  That means we'll try to
                # shrink the widths by a proportion of 192/352--.545454.
                # That would make the first column 48 points, and the second
                # 144 points--adding up to the desired 192.
                #
                # Unfortunately, that's too small for the first column.  It
                # must be 66 points.  Therefore, we go ahead and save that
                # column width as 88 points.  That leaves (192-88==) 104
                # points remaining.  The proportion to shrink the remaining
                # column is (104/264), which, multiplied  by the desired
                # width of 264, is 104: the amount assigned to the remaining
                # column.
                proportion = effectiveRemaining/totalDesired
                # we sort the desired widths by difference between desired and
                # and minimum values, a value called "disappointment" in the
                # code.  This means that the columns with a bigger
                # disappointment will have a better chance of getting more of
                # the available space.
                desiredWidths.sort()
                finalSet = []
                for disappointment, minimum, desired, colNo in desiredWidths:
                    adjusted = proportion * desired
                    if adjusted < minimum:
                        W[colNo] = minimum
                        totalDesired -= desired
                        effectiveRemaining -= minimum
                        if totalDesired:
                            proportion = effectiveRemaining/totalDesired
                    else:
                        finalSet.append((minimum, desired, colNo))
                for minimum, desired, colNo in finalSet:
                    adjusted = proportion * desired
                    assert adjusted >= minimum
                    W[colNo] = adjusted
        else:
            for colNo, minimum in minimums.items():
                W[colNo] = minimum
        if verbose: print 'new widths are:', W
        self._argW = self._colWidths = W
        return W

    def minWidth(self):
        W = list(self._argW)
        width = 0
        elementWidth = self._elementWidth
        rowNos = xrange(self._nrows)
        values = self._cellvalues
        styles = self._cellStyles
        for colNo in xrange(len(W)):
            w = W[colNo]
            if w is None or w=='*' or _endswith(w,'%'):
                final = 0
                for rowNo in rowNos:
                    value = values[rowNo][colNo]
                    style = styles[rowNo][colNo]
                    new = (elementWidth(value, style, (colNo, rowNo)) +
                        style.leftPadding + style.rightPadding)
                    final = max(final, new)
                width += final
            else:
                width += float(w)
        return width # XXX + 1/2*(left and right border widths)

    def _calcSpanRanges(self):
        """
        Work out rects for tables which do row and column spanning.

        This creates some mappings to let the later code determine
        if a cell is part of a "spanned" range.
        self._spanRanges shows the 'coords' in integers of each
        'cell range', or None if it was clobbered:
        (col, row) -> (col0, row0, col1, row1)

        Any cell not in the key is not part of a spanned region.

        This method use absolute data positions so its result can
        be reused after split.
        """
        # Checks if span ranges are already computed.
        if getattr(self, '_spanRanges', None) is not None:
            return

        self._spanRanges = spanRanges = {}
        for x in xrange(self._ncols):
            for y in xrange(self._nrows):
                spanRanges[x,y] = (x, y, x, y)
        self._colSpanCells = []
        self._rowSpanCells = []
        csa = self._colSpanCells.append
        rsa = self._rowSpanCells.append
        for (cmd, start, stop) in self._spanCmds:
            x0, y0 = start
            x1, y1 = stop

            if x0!=x1 or y0!=y1:
                if x0!=x1: #column span
                    for y in xrange(y0, y1+1):
                        for x in xrange(x0,x1+1):
                            csa((x,y))
                if y0!=y1: #row span
                    for y in xrange(y0, y1+1):
                        for x in xrange(x0,x1+1):
                            rsa((x,y))

                for y in xrange(y0, y1+1):
                    for x in xrange(x0,x1+1):
                        spanRanges[x,y] = None
                # set the main entry
                spanRanges[x0,y0] = (x0, y0, x1, y1)

    def _calcNoSplitRanges(self):
        """
        This creates some mappings to let the later code determine
        if a cell is part of a "nosplit" range.
        self._nosplitRanges shows the 'coords' in integers of each
        'cell range', or None if it was clobbered:
        (col, row) -> (col0, row0, col1, row1)

        Any cell not in the key is not part of a spanned region
        """
        # Checks if nosplit ranges are already computed.
        if getattr(self, '_nosplitRanges', None) is not None:
            return

        self._nosplitRanges = nosplitRanges = {}
        for x in xrange(self._ncols):
            for y in xrange(self._nrows):
                nosplitRanges[x,y] = (x, y, x, y)
        self._colNoSplitCells = []
        self._rowNoSplitCells = []
        csa = self._colNoSplitCells.append
        rsa = self._rowNoSplitCells.append
        for (cmd, start, stop) in self._nosplitCmds:
            x0, y0 = start
            x1, y1 = stop

            if x0!=x1 or y0!=y1:
                #column span
                if x0!=x1:
                    for y in xrange(y0, y1+1):
                        for x in xrange(x0,x1+1):
                            csa((x,y))
                #row span
                if y0!=y1:
                    for y in xrange(y0, y1+1):
                        for x in xrange(x0,x1+1):
                            rsa((x,y))

                for y in xrange(y0, y1+1):
                    for x in xrange(x0,x1+1):
                        nosplitRanges[x,y] = None
                # set the main entry
                nosplitRanges[x0,y0] = (x0, y0, x1, y1)

    def _calcSpanRects(self):
        """
        Work out rects for tables which do row and column spanning.

        Based on self._spanRanges, which is already known,
        and the widths which were given or previously calculated,
        self._spanRects shows the real coords for drawing:

            (col, row) -> (x, y, width, height)

        for each cell.  Any cell which 'does not exist' as another
        has spanned over it will get a None entry on the right.

        This method generates relative positions so its results cannot
        be reused after split.
        """
        if getattr(self,'_spanRects',None): return
        colpositions = self._colpositions
        rowpositions = self._rowpositions
        self._spanRects = spanRects = {}
        self._vBlocks = vBlocks = {}
        self._hBlocks = hBlocks = {}
        for (coord, value) in self._spanRanges.items():
            if value is None:
                spanRects[coord] = None
            else:
                col,row = coord
                # Testing row for visibility should be enough since no splits
                # are permitted across spanned areas.
                if not self._is_visible_row(row):
                    continue
                col0, row0, col1, row1 = value
                row0 = self._abs_to_vis(row0)
                row1 = self._abs_to_vis(row1)
                if col1-col0>0:
                    for _ in xrange(col0+1,col1+1):
                        vBlocks.setdefault(colpositions[_],[]).append((rowpositions[row1+1],rowpositions[row0]))
                if row1-row0>0:
                    for _ in xrange(row0+1,row1+1):
                        hBlocks.setdefault(rowpositions[_],[]).append((colpositions[col0],colpositions[col1+1]))
                x = colpositions[col0]
                y = rowpositions[row1+1]
                width = colpositions[col1+1] - x
                height = rowpositions[row0] - y
                spanRects[coord] = (x, y, width, height)

        for _ in hBlocks, vBlocks:
            for value in _.values():
                value.sort()

    def setStyle(self, tblstyle):
        if not isinstance(tblstyle,TableStyle):
            tblstyle = TableStyle(tblstyle)
        for cmd in tblstyle.getCommands():
            self._addCommand(cmd)
        for k,v in tblstyle._opts.items():
            setattr(self,k,v)
        for a in ('spaceBefore','spaceAfter'):
            if not hasattr(self,a) and hasattr(tblstyle,a):
                setattr(self,a,getattr(tblstyle,a))

    def _normalizeCoord(self, sc, ec, sr, er):
        """
        Normalizes cols/rows coordinates.
        """
        if sc < 0: sc = sc + self._ncols
        if ec < 0: ec = ec + self._ncols
        if sr < 0: sr = sr + self._nrows
        if er < 0: er = er + self._nrows
        return sc, ec, sr, er

    def _addCommand(self,cmd):
        if cmd[0] in ('BACKGROUND','ROWBACKGROUNDS','COLBACKGROUNDS'):
            op, (sc, sr), (ec, er), arg = cmd
            sc, ec, sr, er = self._normalizeCoord(sc, ec, sr, er)
            cmd = (op, (sc, sr), (ec, er), arg)
            self._bkgrndcmds.append(cmd)
        elif cmd[0] == 'SPAN':
            op, (sc, sr), (ec, er) = cmd
            sc, ec, sr, er = self._normalizeCoord(sc, ec, sr, er)
            if sc > ec: sc, ec = ec, sc
            if sr > er: sr, er = er, sr
            cmd = (op, (sc, sr), (ec, er))
            self._spanCmds.append(cmd)
        elif cmd[0] == 'NOSPLIT':
            op, (sc, sr), (ec, er) = cmd
            sc, ec, sr, er = self._normalizeCoord(sc, ec, sr, er)
            if sc > ec: sc, ec = ec, sc
            if sr > er: sr, er = er, sr
            cmd = (op, (sc, sr), (ec, er))
            self._nosplitCmds.append(cmd)
        elif _isLineCommand(cmd):
            # we expect op, start, stop, weight, colour, cap, dashes, join
            cmd = list(cmd)
            if len(cmd)<5: raise ValueError('bad line command '+str(cmd))

            #determine line cap value at position 5. This can be str or numeric.
            if len(cmd)<6:
                cmd.append(1)
            else:
                cap = _convert2int(cmd[5], LINECAPS, 0, 2, 'cap', cmd)
                cmd[5] = cap

            #dashes at index 6 - this is a dash array:
            if len(cmd)<7: cmd.append(None)

            #join mode at index 7 - can be str or numeric, look up as for caps
            if len(cmd)<8: cmd.append(1)
            else:
                join = _convert2int(cmd[7], LINEJOINS, 0, 2, 'join', cmd)
                cmd[7] = join

            #linecount at index 8.  Default is 1, set to 2 for double line.
            if len(cmd)<9: cmd.append(1)
            else:
                lineCount = cmd[8]
                if lineCount is None:
                    lineCount = 1
                    cmd[8] = lineCount
                assert lineCount >= 1
            #linespacing at index 9. Not applicable unless 2+ lines, defaults to line
            #width so you get a visible gap between centres
            if len(cmd)<10: cmd.append(cmd[3])
            else:
                space = cmd[9]
                if space is None:
                    space = cmd[3]
                    cmd[9] = space
            assert len(cmd) == 10

            (op, (sc,sr), (ec,er), weight, color, cap, dash, join, count,
                space) = cmd[:]
            sc, ec, sr, er = self._normalizeCoord(sc, ec, sr, er)
            cmd = (op, (sc,sr), (ec,er), weight, color, cap, dash, join, count,
                space)
            self._linecmds.append(cmd)
        else:
            (op, (sc, sr), (ec, er)), values = cmd[:3] , cmd[3:]
            sc, ec, sr, er = self._normalizeCoord(sc, ec, sr, er)
            for i in xrange(sr, er+1):
                for j in xrange(sc, ec+1):
                    _setCellStyle(self._cellStyles, i, j, op, values)

    def _drawLines(self):
        ccap, cdash, cjoin = None, None, None
        self.canv.saveState()
        for op, (sc,sr), (ec,er), weight, color, cap, dash, join, count, space in self._linecmds:
            if isinstance(sr,basestring) and sr.startswith('split'): continue
            if cap!=None and ccap!=cap:
                self.canv.setLineCap(cap)
                ccap = cap
            if dash is None or dash == []:
                if cdash is not None:
                    self.canv.setDash()
                    cdash = None
            elif dash != cdash:
                self.canv.setDash(dash)
                cdash = dash
            if join is not None and cjoin!=join:
                self.canv.setLineJoin(join)
                cjoin = join
            getattr(self,_LineOpMap.get(op, '_drawUnknown' ))( (sc, sr), (ec, er), weight, color, count, space)
        self.canv.restoreState()
        self._curcolor = None

    def _drawUnknown(self,  (sc, sr), (ec, er), weight, color, count, space):
        #we are only called from _drawLines which is one level up
        import sys
        op = sys._getframe(1).f_locals['op']
        raise ValueError("Unknown line command '%s'" % op)

    def _is_visible_line(self, line_num):
        """
        Checks if line is in visible area of table.
        """
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        return (line_num <= self.repeatRows or
            activeRows0 <= line_num <= self._activeRows[1] or
            self._nrows - self._repeatRowsB <= line_num <= self._nrows)

    def _is_visible_row(self, row_num):
        """
        Checks if row is in visible area of table.
        """
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        return (row_num < self.repeatRows or
            activeRows0 <= row_num < self._activeRows[1] or
            self._nrows - self._repeatRowsB <= row_num < self._nrows)

    def _abs_to_vis(self, line_num):
        """
        Translates absolute line positions to relative (visible positions).
        """
        if line_num <= self.repeatRows:
            return line_num
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        if activeRows0 <= line_num <= self._activeRows[1]:
            line_num -= activeRows0 - self.repeatRows
            return line_num
        if self._nrows - self._repeatRowsB <= line_num <= self._nrows:
            line_num -= activeRows0 - self.repeatRows
            line_num -= self._nrows - self._repeatRowsB - self._activeRows[1]
            return line_num
        raise ValueError('`line_num` outside visible area!')

    def _vis_to_abs(self, line_num):
        """
        Translates relative line positions to absolute.
        """
        if line_num <= self.repeatRows:
            return line_num
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        line_num += activeRows0 - self.repeatRows
        if activeRows0 <= line_num <= self._activeRows[1]:
            return line_num
        line_num += self._nrows - self._repeatRowsB - self._activeRows[1]
        if self._nrows - self._repeatRowsB <= line_num <= self._nrows:
            return line_num
        raise ValueError('`line_num` outside visible area!')

    def _drawGrid(self, (sc, sr), (ec, er), weight, color, count, space):
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        # Checks if whole grid is outside visible area.
        if sr >= self.repeatRows and er < activeRows0:
            return
        if sr >= self._activeRows[1] and er < self._nrows - self._repeatRowsB:
            return

        # Some parts visible - rendering.
        self._drawBox( (sc, sr), (ec, er), weight, color, count, space)
        self._drawInnerGrid( (sc, sr), (ec, er), weight, color, count, space)

    def _drawBox(self,  (sc, sr), (ec, er), weight, color, count, space):
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        # Checks if whole box is outside visible area.
        if sr >= self.repeatRows and er < activeRows0:
            return
        if sr >= self._activeRows[1] and er < self._nrows - self._repeatRowsB:
            return

        # Some parts visible - rendering.

        # If start row visible renders upper horizontal line.
        if self._is_visible_row(sr):
            self._drawHLines((sc, sr), (ec, sr), weight, color, count, space)
        # If end row visible renders lower horizontal line.
        if self._is_visible_row(er):
            self._drawHLines((sc, er+1), (ec, er+1), weight, color, count, space)

        # Renders vertical lines.
        self._drawVLines((sc, sr), (sc, er + 1), weight, color, count, space)
        self._drawVLines((ec+1, sr), (ec+1, er + 1), weight, color, count, space)

    def _drawInnerGrid(self, (sc, sr), (ec, er), weight, color, count, space):
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        # Checks if whole inner grid is outside visible area.
        if sr >= self.repeatRows and er < activeRows0:
            return
        if sr >= self._activeRows[1] and er < self._nrows - self._repeatRowsB:
            return

        # Some parts visible - rendering.
        self._drawHLines((sc, sr+1), (ec, er), weight, color, count, space)
        self._drawVLines((sc+1, sr), (ec, er + 1), weight, color, count, space)

    def _drawLineAbove(self, (sc, sr), (ec, er), weight, color, count, space):
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        # Checks if whole row range is outside visible area.
        if sr >= self.repeatRows and er < activeRows0:
            return
        if sr >= self._activeRows[1] and er < self._nrows - self._repeatRowsB:
            return

        # Some parts visible - searching for visible rows.
        visible = []
        for i in xrange(sr, er + 1):
            if not self._is_visible_row(i):
                continue
            visible.append(i)

        # Generates line for each visible row.
        for vis in visible:
            self._drawHLines((sc, vis), (ec, vis), weight, color, count, space)

    def _drawLineBelow(self, (sc, sr), (ec, er), weight, color, count, space):
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        # Checks if whole row range is outside visible area.
        if sr >= self.repeatRows and er < activeRows0:
            return
        if sr >= self._activeRows[1] and er < self._nrows - self._repeatRowsB:
            return

        # Some parts visible - searching for visible rows.
        visible = []
        for i in xrange(sr, er + 1):
            if not self._is_visible_row(i):
                continue
            visible.append(i)

        # Generates line for each visible row.
        for vis in visible:
            self._drawHLines((sc, vis + 1), (ec, vis +1), weight, color, count, space)

    def _prepLine(self, weight, color):
        if color != self._curcolor:
            self.canv.setStrokeColor(color)
            self._curcolor = color
        if weight != self._curweight:
            self.canv.setLineWidth(weight)
            self._curweight = weight

    def _drawHLines(self, (sc, sr), (ec, er), weight, color, count, space):
        ecp = self._colpositions[sc:ec+2]

        visible = [i for i in xrange(sr, er + 1) if self._is_visible_line(i)]
        rp_pos = [self._abs_to_vis(abs_num) for abs_num in visible]
        rp = [self._rowpositions[pos] for pos in rp_pos]

        if len(ecp)<=1 or len(rp)<1: return
        self._prepLine(weight, color)
        scp = ecp[0]
        ecp = ecp[-1]
        hBlocks = getattr(self,'_hBlocks',{})
        canvLine = self.canv.line
        if count == 1:
            for y in rp:
                _hLine(canvLine, scp, ecp, y, hBlocks)
        else:
            lf = lambda x0,y0,x1,y1,canvLine=canvLine, ws=weight+space, count=count: _multiLine(x0,x1,y0,canvLine,ws,count)
            for y in rp:
                _hLine(lf, scp, ecp, y, hBlocks)

    def _drawVLines(self, (sc, sr), (ec, er), weight, color, count, space):
        visible = [i for i in xrange(sr, er + 1) if self._is_visible_line(i)]
        rp_pos = [self._abs_to_vis(abs_num) for abs_num in visible]
        erp = [self._rowpositions[pos] for pos in rp_pos]

        cp  = self._colpositions[sc:ec+1]
        if len(erp)<=1 or len(cp)<1: return
        self._prepLine(weight, color)
        srp = erp[0]
        erp = erp[-1]
        vBlocks = getattr(self,'_vBlocks',{})
        canvLine = lambda y0, x0, y1, x1, _line=self.canv.line: _line(x0,y0,x1,y1)
        if count == 1:
            for x in cp:
                _hLine(canvLine, erp, srp, x, vBlocks)
        else:
            lf = lambda x0,y0,x1,y1,canvLine=canvLine, ws=weight+space, count=count: _multiLine(x0,x1,y0,canvLine,ws,count)
            for x in cp:
                _hLine(lf, erp, srp, x, vBlocks)

    def _drawLineAfter(self, (sc, sr), (ec, er), weight, color, count, space):
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        # Checks if whole row range is outside visible area.
        if sr >= self.repeatRows and er + 1 < activeRows0:
            return
        if sr >= self._activeRows[1] and er + 1 < self._nrows - self._repeatRowsB:
            return

        # Some parts visible - rendering.
        self._drawVLines((sc + 1, sr), (ec + 1, er + 1), weight, color, count, space)

    def _drawLineBefore(self, (sc, sr), (ec, er), weight, color, count, space):
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        # Checks if whole row range is outside visible area.
        if sr >= self.repeatRows and er + 1 < activeRows0:
            return
        if sr >= self._activeRows[1] and er + 1 < self._nrows - self._repeatRowsB:
            return

        # Some parts visible - rendering.
        self._drawVLines((sc, sr), (ec, er + 1), weight, color, count, space)

    def wrap(self, availWidth, availHeight):
        self._calc(availWidth, availHeight)
        #nice and easy, since they are predetermined size
        self.availWidth = availWidth
        return (self._width, self._height)

    def onSplit(self,T,byRow=1):
        '''
        This method will be called when the Table is split.
        Special purpose tables can override to do special stuff.
        '''
        pass

    def _splitRows(self,availHeight):
        n=self._getFirstPossibleSplitRowPosition(availHeight)
        if n==0: return []
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        lim = len(self._rowHeights[activeRows0:self._activeRows[1]])
        if n==lim: return [self]

        R0 = self._copy()
        # We set what part of data is going to be visible.
        R0._activeRows = (activeRows0, activeRows0 + n)

        R1 = self._copy()
        # We set what part of data is going to be visible.
        R1._activeRows = (activeRows0 + n, self._activeRows[1])

        self.onSplit(R0)
        self.onSplit(R1)
        return [R0,R1]

    def _getRowImpossible(impossible,cells,ranges):
        """
        Marks row numbers where split is impossible due to span or nosplit
        commands.
        """
        for xy in cells:
            r=ranges[xy]
            if r!=None:
                y1,y2=r[1],r[3]
                if y1!=y2:
                    ymin=min(y1,y2) #normalize
                    ymax=max(y1,y2) #normalize
                    y=ymin
                    while 1:
                        if y>=ymax: break
                        impossible[y]=None #split at position y is impossible because of overlapping rowspan
                        y+=1
    _getRowImpossible=staticmethod(_getRowImpossible)

    def _getFirstPossibleSplitRowPosition(self,availHeight):
        # Since impossible do not change after split, we try to reuse old data
        impossible = getattr(self, '_impossible', None)
        if impossible is None:
            impossible={}
            if self._spanCmds:
                self._getRowImpossible(impossible, self._rowSpanCells,
                    self._spanRanges)
            if self._nosplitCmds:
                self._getRowImpossible(impossible, self._rowNoSplitCells,
                    self._nosplitRanges)
            self._impossible = impossible

        h = sum(self._rowHeights[:self.repeatRows])
        h += sum(self._rowHeights[self._nrows-self._repeatRowsB:])
        split_at = 0 # from this point of view 0 is the first position where the table may *always* be splitted
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        for n, rh in enumerate(self._rowHeights[activeRows0:self._activeRows[1]]):
            if h+rh>availHeight:
                break
            if not impossible.has_key(self._vis_to_abs(n+self.repeatRows)):
                split_at=n + 1
            h=h+rh
        return split_at

    def split(self, availWidth, availHeight):
        self._calc(availWidth, availHeight)
        if self.splitByRow:
            if not rl_config.allowTableBoundsErrors and self._width>availWidth: return []
            return self._splitRows(availHeight)
        else:
            raise NotImplementedError

    def draw(self):
        self._curweight = self._curcolor = self._curcellstyle = None
        self._drawBkgrnd()
        activeRows0 = self._activeRows[0] if self._activeRows[0] is not None else self.repeatRows # ugly hack to make it backward compatible
        row_nums = (range(0, self.repeatRows) +
            range(activeRows0, self._activeRows[1]) +
            range(self._nrows-self._repeatRowsB, self._nrows))
        if not self._spanCmds:
            # old fashioned case, no spanning, steam on and do each cell
            cellvalues = (self._cellvalues[:self.repeatRows] +
                self._cellvalues[activeRows0:self._activeRows[1]] +
                self._cellvalues[self._nrows-self._repeatRowsB:])
            cellstyles = (self._cellStyles[:self.repeatRows] +
                self._cellStyles[activeRows0:self._activeRows[1]] +
                self._cellStyles[self._nrows-self._repeatRowsB:])
            rowheights = (self._rowHeights[:self.repeatRows] +
                self._rowHeights[activeRows0:self._activeRows[1]] +
                self._rowHeights[self._nrows-self._repeatRowsB:])
            for rowNo, row, rowstyle, rowpos, rowheight in zip(row_nums,
                cellvalues, cellstyles, self._rowpositions[1:], rowheights):

                for colNo, cellval, cellstyle, colpos, colwidth in zip(
                    xrange(self._ncols), row, rowstyle,
                    self._colpositions[:-1], self._colWidths):

                    if isinstance(cellval, Formula):
                        cellval = cellval(self._cellvalues, self.repeatRows,
                            self._repeatRowsB,
                            (activeRows0, self._activeRows[1]),
                            (colNo, rowNo))
                    self._drawCell(cellval, cellstyle, (colpos, rowpos),
                        (colwidth, rowheight))
        else:
            # we have some row or col spans, need a more complex algorithm
            # to find the rect for each
            for rowNo in row_nums:
                for colNo in xrange(self._ncols):
                    cellRect = self._spanRects[colNo, rowNo]
                    if cellRect is not None:
                        (x, y, width, height) = cellRect
                        cellval = self._cellvalues[rowNo][colNo]
                        cellstyle = self._cellStyles[rowNo][colNo]
                        if isinstance(cellval, Formula):
                            cellval = cellval(self._cellvalues,
                                self.repeatRows, self._repeatRowsB,
                                (activeRows0, self._activeRows[1]),
                                (colNo, rowNo))
                        self._drawCell(cellval, cellstyle, (x, y),
                            (width, height))
        self._drawLines()

    def _drawBkgrnd(self):
        nrows = self._nrows
        ncols = self._ncols
        canv = self.canv
        colpositions = self._colpositions
        rowpositions = self._rowpositions
        rowHeights = self._rowHeights
        colWidths = self._colWidths
        spanRects = getattr(self,'_spanRects',None)
        for cmd, (sc, sr), (ec, er), arg in self._bkgrndcmds:

            visible = []
            for row_num in xrange(sr, er + 1):
                if not self._is_visible_row(row_num):
                    continue
                visible.append(row_num)
            if not visible:
                continue
            sr = self._abs_to_vis(visible[0])
            er = self._abs_to_vis(visible[-1])

            x0 = colpositions[sc]
            y0 = rowpositions[sr]
            x1 = colpositions[min(ec+1,ncols)]
            y1 = rowpositions[er+1]
            w, h = x1-x0, y1-y0
            if callable(arg):
                arg(self,canv, x0, y0, w, h)
            elif cmd == 'ROWBACKGROUNDS':
                #Need a list of colors to cycle through.  The arguments
                #might be already colours, or convertible to colors, or
                # None, or the str 'None'.
                #It's very common to alternate a pale shade with None.
                colorCycle = map(colors.toColorOrNone, arg)
                count = len(colorCycle)
                for i, abs_row_num in enumerate(visible):
                    color = colorCycle[i%count]
                    h = rowHeights[abs_row_num]
                    if color:
                        canv.setFillColor(color)
                        canv.rect(x0, y0, w, -h, stroke=0,fill=1)
                    y0 = y0 - h
            elif cmd == 'COLBACKGROUNDS':
                #cycle through colours columnwise
                colorCycle = map(colors.toColorOrNone, arg)
                count = len(colorCycle)
                colCount = ec - sc + 1
                for i in xrange(colCount):
                    color = colorCycle[i%count]
                    w = colWidths[sc + i]
                    if color:
                        canv.setFillColor(color)
                        canv.rect(x0, y0, w, h, stroke=0,fill=1)
                    x0 = x0 +w
            else:   #cmd=='BACKGROUND'
                color = colors.toColorOrNone(arg)
                if color:
                    if ec==sc and er==sr and spanRects:
                        xywh = spanRects.get((sc,sr))
                        if xywh:
                            #it's a single cell
                            x0, y0, w, h = xywh
                    canv.setFillColor(color)
                    canv.rect(x0, y0, w, h, stroke=0,fill=1)

    def _drawCell(self, cellval, cellstyle, (colpos, rowpos), (colwidth, rowheight)):
        if self._curcellstyle is not cellstyle:
            cur = self._curcellstyle
            if cur is None or cellstyle.color != cur.color:
                self.canv.setFillColor(cellstyle.color)
            if cur is None or cellstyle.leading != cur.leading or cellstyle.fontname != cur.fontname or cellstyle.fontsize != cur.fontsize:
                self.canv.setFont(cellstyle.fontname, cellstyle.fontsize, cellstyle.leading)
            self._curcellstyle = cellstyle

        just = cellstyle.alignment
        valign = cellstyle.valign
        if isinstance(cellval,(tuple,list,Flowable)):
            if not isinstance(cellval,(tuple,list)): cellval = (cellval,)
            # we assume it's a list of Flowables
            W = []
            H = []
            w, h = self._listCellGeom(cellval,colwidth,cellstyle,W=W, H=H,aH=rowheight)
            if valign=='TOP':
                y = rowpos + rowheight - cellstyle.topPadding
            elif valign=='BOTTOM':
                y = rowpos+cellstyle.bottomPadding + h
            else:
                y = rowpos+(rowheight+cellstyle.bottomPadding-cellstyle.topPadding+h)/2.0
            if cellval: y += cellval[0].getSpaceBefore()
            for v, w, h in zip(cellval,W,H):
                if just=='LEFT': x = colpos+cellstyle.leftPadding
                elif just=='RIGHT': x = colpos+colwidth-cellstyle.rightPadding - w
                elif just in ('CENTRE', 'CENTER'):
                    x = colpos+(colwidth+cellstyle.leftPadding-cellstyle.rightPadding-w)/2.0
                else:
                    raise ValueError('Invalid justification %s' % just)
                y -= v.getSpaceBefore()
                y -= h
                v.drawOn(self.canv,x,y)
                y -= v.getSpaceAfter()
        else:
            if just == 'LEFT':
                draw = self.canv.drawString
                x = colpos + cellstyle.leftPadding
            elif just in ('CENTRE', 'CENTER'):
                draw = self.canv.drawCentredString
                x = colpos+(colwidth+cellstyle.leftPadding-cellstyle.rightPadding)*0.5
            elif just == 'RIGHT':
                draw = self.canv.drawRightString
                x = colpos + colwidth - cellstyle.rightPadding
            elif just == 'DECIMAL':
                draw = self.canv.drawAlignedString
                x = colpos + colwidth - cellstyle.rightPadding
            else:
                raise ValueError('Invalid justification %s' % just)
            vals = str(cellval).split("\n")
            n = len(vals)
            leading = cellstyle.leading
            fontsize = cellstyle.fontsize
            if valign=='BOTTOM':
                y = rowpos + cellstyle.bottomPadding+n*leading-fontsize
            elif valign=='TOP':
                y = rowpos + rowheight - cellstyle.topPadding - fontsize
            elif valign=='MIDDLE':
                #tim roberts pointed out missing fontsize correction 2004-10-04
                y = rowpos + (cellstyle.bottomPadding + rowheight-cellstyle.topPadding+n*leading)/2.0 - fontsize
            else:
                raise ValueError("Bad valign: '%s'" % str(valign))

            for v in vals:
                draw(x, y, v)
                y -= leading
            onDraw = getattr(cellval,'onDraw',None)
            if onDraw:
                onDraw(self.canv,cellval.kind,cellval.label)

        if cellstyle.href:
            #external hyperlink
            self.canv.linkURL(cellstyle.href, (colpos, rowpos, colpos + colwidth, rowpos + rowheight), relative=1)
        if cellstyle.destination:
            #external hyperlink
            self.canv.linkRect("", cellstyle.destination, Rect=(colpos, rowpos, colpos + colwidth, rowpos + rowheight), relative=1)

    def _copy(self):
        """
        Makes copy of self.
        """
        shadow = self.__class__(data = self._cellvalues,
            colWidths = self._colWidths, rowHeights = self._rowHeights,
            repeatRows = self.repeatRows, repeatCols = self.repeatCols,
            splitByRow = self.splitByRow, normalizedData = 1,
            cellStyles = self._cellStyles, activeRows = self._activeRows,
            repeatRowsB = self._repeatRowsB, hAlign = self.hAlign,
            vAlign = self.vAlign, ident = self.ident)

        # copy the commands
        shadow._linecmds = self._linecmds
        shadow._bkgrndcmds = self._bkgrndcmds
        shadow._spanCmds = self._spanCmds
        shadow._nosplitCmds = self._nosplitCmds

        # copy span related data
        if getattr(self, '_spanRanges', None) is not None:
            shadow._spanRanges = self._spanRanges
            shadow._colSpanCells = self._colSpanCells
            shadow._rowSpanCells = self._rowSpanCells

        # copy nosplit related data
        if getattr(self, '_nosplitRanges', None) is not None:
            shadow._nosplitRanges = self._nosplitRanges
            shadow._colNoSplitCells = self._colNoSplitCells
            shadow._rowNoSplitCells = self._rowNoSplitCells

        shadow._impossible = getattr(self, '_impossible', None)

        return shadow


_LineOpMap = {  'GRID':'_drawGrid',
                'BOX':'_drawBox',
                'OUTLINE':'_drawBox',
                'INNERGRID':'_drawInnerGrid',
                'LINEBELOW':'_drawLineBelow',
                'LINEABOVE':'_drawLineAbove',
                'LINEBEFORE':'_drawLineBefore',
                'LINEAFTER':'_drawLineAfter', }


if __name__ == '__main__':
    from demo import demo

    demo()
