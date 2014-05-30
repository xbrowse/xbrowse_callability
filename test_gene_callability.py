from nose.tools import assert_equals
from nose.tools import assert_raises
from nose.tools import assert_items_equal
import gene_callability as gc


#######################################################################################################################
# Tests for _get_missing_interval_list
#######################################################################################################################
def test_get_features_with_missing_intervals():
    features = [{'start': 1000, 'end': 2000}, {'start': 2000, 'end': 3000}, {'start': 2500, 'end': 4000}]
    qmi_results = [{'INTERVAL': '1:50-100'}, {'INTERVAL': '1:700-1050'}]
    result = [({'start': 1000, 'end': 2000}, [{'INTERVAL': '1:700-1050'}])]

    assert_items_equal(result, gc.get_features_with_missing_intervals(features, qmi_results))


#######################################################################################################################
# Tests for _parse_interval
#######################################################################################################################
def test_parse_interval_with_standard_input():
    input = "1:34636-40000"
    result = ('1', 34636, 40000)

    assert_equals(result, gc._parse_interval(input))


def test_parse_interval_with_single_base_interval():
    input = "2:34636"
    result = ('2', 34636, 34636)

    assert_equals(result, gc._parse_interval(input))


def test_parse_interval_with_invalid_format():
    input = "1>34636+40000"

    assert_raises(ValueError, gc._parse_interval, input)


#######################################################################################################################
# Tests for _IntervalOverlapCalculator
#######################################################################################################################
class TestIntervalOverlapCalculator():

    def setup(self):
        self.interval_tree = gc._IntervalOverlapCalculator([(1, 5), (10, 20), (15, 25), (40, 50)])

    def test_overlap_with_one_overlapping_interval(self):
        input = (2, 6)
        result = [(1, 5)]
        assert_items_equal(result, self.interval_tree.get_overlapping_intervals(input))

    def test_overlap_with_multiple_overlapping_intervals(self):
        input = (10, 25)
        result = [(10, 20), (15, 25)]
        assert_items_equal(result, self.interval_tree.get_overlapping_intervals(input))


#######################################################################################################################
# Tests for _get_missing_interval_map
#######################################################################################################################
def test_get_missing_interval_map():
    input = [{'INTERVAL': '1:1000-2000', 'other': 'other_data'}, {'INTERVAL': '1:2000-3000', 'other': 'other_data'}]
    result = {1000: {'INTERVAL': '1:1000-2000', 'other':'other_data'}, 2000: {'INTERVAL': '1:2000-3000', 'other':'other_data'}}

    assert_equals(result, gc._get_missing_interval_map(input))


#######################################################################################################################
# Tests for _get_missing_interval_list
#######################################################################################################################
def test_get_missing_interval_list():
    input = [{'INTERVAL': '1:1000-2000', 'other': 'other_data'}, {'INTERVAL': '1:2000-3000', 'other': 'other_data'}]
    result = [(1000, 2000), (2000, 3000)]

    assert_equals(result, gc._get_missing_interval_list(input))