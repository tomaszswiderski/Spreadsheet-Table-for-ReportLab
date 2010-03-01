# -*- coding: utf-8 -*-

__author__ = u'Tomasz Świderski <contact@tomaszswiderski.com>'
__copyright__ = u'Copyright (c) 2010 Tomasz Świderski'

from decimal import Decimal as D, InvalidOperation as DConversionError


class Formula(object):

    """
    Abstract class of all formulas.
    """

    def __init__(self, longest_value=None):
        # This value will be used to determine cell width during Table's
        # self._calc_width.
        self._longest_value = longest_value

    def get_max_value(self, data, repeat_rows, repeat_rows_b, cell_coord):
        """
        Returns largest possible value based on data and coord information.
        """
        if self._longest_value is not None:
            return self._longest_value
        return self._get_max_value(data, repeat_rows, repeat_rows_b,
            cell_coord)

    def __call__(self, data, repeat_rows, repeat_rows_b, active_rows,
            cell_coord):
        """
        Return evaluated string value.
        """
        raise NotImplementedError

    def _get_max_value(self, data, repeat_rows, repeat_rows_b, cell_coord):
        """
        Returns largest possible value of Formula.
        """
        raise NotImplementedError


class CurrentPageColSum(Formula):

    """
    Calculates sum of column from current page only.
    """

    def __init__(self, decimal_places=2, longest_value=None,
            ignore_convert_errors=True):
        Formula.__init__(self, longest_value)
        self._ignore_convert_errors = ignore_convert_errors
        self._decimal_places = decimal_places

    def __call__(self, data, repeat_rows, repeat_rows_b, active_rows,
            cell_coord):
        """
        Calculates sum of column. Only values inside current active rows will
        be considered.  All data will be converted to decimals. Unconvertable
        values can be ignored or cause to raise Exception.
        """
        # Takes slice of data which will be used in evaluation.
        cell_row = cell_coord[1]
        if active_rows[0] <= cell_row < active_rows[1]:
            raise ValueError('Formula inside range to be evaluated!')

        active_data = data[active_rows[0]:active_rows[1]]

        col_num = cell_coord[0]
        sum = D(0)

        for row_num, row in enumerate(active_data):
            col_value = row[col_num]
            if isinstance(col_value, Formula):
                col_value = col_value(data, repeat_rows, repeat_rows_b,
                    active_rows, (col_num, row_num + active_rows[0]))
            try:
                value = D(col_value)
            except DConversionError:
                if not self._ignore_convert_errors:
                    raise
                continue
            sum += value
        format_str = '%%.%df' % self._decimal_places

        return format_str % sum

    def _get_max_value(self, data, repeat_rows, repeat_rows_b, cell_coord):
        """
        Returns largest possible value of Formula.
        """
        # Takes slice of data which will be used in evaluation.
        cell_row = cell_coord[1]
        end_row = len(data) - repeat_rows_b
        if repeat_rows <= cell_row < end_row:
            raise ValueError('Formula inside range to be evaluated!')

        active_data = data[repeat_rows:end_row]

        col_num = cell_coord[0]
        sum = D(0)

        for row_num, row in enumerate(active_data):
            col_value = row[col_num]
            if isinstance(col_value, Formula):
                col_value = col_value(data, repeat_rows, repeat_rows_b,
                    (repeat_rows, end_row), (col_num, row_num + repeat_rows))
            try:
                value = D(col_value)
            except DConversionError:
                if not self._ignore_convert_errors:
                    raise
                continue
            sum += value
        format_str = '%%.%df' % self._decimal_places

        return format_str % sum


class PreviousPagesColSum(Formula):

    """
    Calculates sum of column from previous pages.
    """

    def __init__(self, starting_value=D(0), decimal_places=2,
            ignore_convert_errors=True, longest_value=None):
        Formula.__init__(self, longest_value)
        self._ignore_convert_errors = ignore_convert_errors
        self._starting_value = starting_value
        self._decimal_places = decimal_places

    def __call__(self, data, repeat_rows, repeat_rows_b, active_rows,
            cell_coord):
        """
        Calculates sum of column from beginning to start of current page.
        All data will be converted to decimals. Unconvertable values
        can be ignored or cause to raise Exception.
        """
        # Takes slice of data which will be used in evaluation.
        cell_row = cell_coord[1]
        if repeat_rows <= cell_row < active_rows[0]:
            raise ValueError('Formula inside range to be evaluated!')

        active_data = data[repeat_rows:active_rows[0]]

        col_num = cell_coord[0]
        sum = self._starting_value

        for row_num, row in enumerate(active_data):
            col_value = row[col_num]
            if isinstance(col_value, Formula):
                col_value = col_value(data, repeat_rows, repeat_rows_b,
                    active_rows, (col_num, row_num + active_rows[0]))
            try:
                value = D(col_value)
            except DConversionError:
                if not self._ignore_convert_errors:
                    raise
                continue
            sum += value
        format_str = '%%.%df' % self._decimal_places

        return format_str % sum

    def _get_max_value(self, data, repeat_rows, repeat_rows_b, cell_coord):
        """
        Returns largest possible value of Formula.
        """
        # Takes slice of data which will be used in evaluation.
        cell_row = cell_coord[1]
        end_row = len(data) - repeat_rows_b
        if repeat_rows <= cell_row < end_row:
            raise ValueError('Formula inside range to be evaluated!')

        active_data = data[repeat_rows:end_row]

        col_num = cell_coord[0]
        sum = self._starting_value

        for row_num, row in enumerate(active_data):
            col_value = row[col_num]
            if isinstance(col_value, Formula):
                col_value = col_value(data, repeat_rows, repeat_rows_b,
                    (repeat_rows, end_row), (col_num, row_num + repeat_rows))
            try:
                value = D(col_value)
            except DConversionError:
                if not self._ignore_convert_errors:
                    raise
                continue
            sum += value
        format_str = '%%.%df' % self._decimal_places

        return format_str % sum


class TotalPagesColSum(Formula):

    """
    Calculates sum of column from previous pages and current page.
    """

    def __init__(self, starting_value=D(0), decimal_places=2,
            ignore_convert_errors=True, longest_value=None):
        Formula.__init__(self, longest_value)
        self._prev_pages_col_sum = PreviousPagesColSum(
            starting_value = starting_value, decimal_places = decimal_places,
            ignore_convert_errors = ignore_convert_errors,
            longest_value = longest_value)
        self._curr_page_col_sum = CurrentPageColSum(
            decimal_places = decimal_places,
            ignore_convert_errors = ignore_convert_errors,
            longest_value = longest_value)
        self._decimal_places = decimal_places

    def __call__(self, data, repeat_rows, repeat_rows_b, active_rows,
            cell_coord):
        """
        Calculates sum of column up to current formula coordinates.
        All data will be converted to decimals. Unconvertable values
        can be ignored or cause to raise Exception.
        """
        prev_value = self._prev_pages_col_sum(data, repeat_rows, repeat_rows_b,
            active_rows, cell_coord)
        curr_value = self._curr_page_col_sum(data, repeat_rows, repeat_rows_b,
            active_rows, cell_coord)
        sum = D(prev_value) + D(curr_value)
        format_str = '%%.%df' % self._decimal_places

        return format_str % sum

    def _get_max_value(self, data, repeat_rows, repeat_rows_b, cell_coord):
        """
        Returns largest possible value of Formula.
        """
        return self._prev_pages_col_sum.get_max_value(data, repeat_rows,
            repeat_rows_b, cell_coord)


class RowNumber(Formula):

    """
    Returns number of current row starting from beginning of current page.
    """

    def __init__(self, starting_value=1, longest_value=None):
        Formula.__init__(self, longest_value)
        self._starting_value = starting_value

    def __call__(self, data, repeat_rows, repeat_rows_b, active_rows,
            cell_coord):
        """
        Returns row number of formula.
        """
        cell_row = cell_coord[1]
        if not active_rows[0] <= cell_row < active_rows[1]:
            raise ValueError('Formula must be inside visible range!')

        return str(cell_row - active_rows[0] + self._starting_value)

    def _get_max_value(self, data, repeat_rows, repeat_rows_b, cell_coord):
        """
        Returns largest possible value of Formula.
        """
        # Takes slice of data which will be used in evaluation.
        cell_row = cell_coord[1]
        end_row = len(data) - repeat_rows_b
        if not repeat_rows <= cell_row < end_row:
            raise ValueError('Formula must be inside visible range!')

        active_data = data[repeat_rows:end_row]

        return str(len(active_data))
