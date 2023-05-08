import pytest
from asimtools.utils import join_names

# @pytest.fixture
# def input_value():
#    input = 39
#    return input

@pytest.mark.parametrize("input, expected",[
    (['a', 'b'],'a_b'),
    (['a', 'b', 'c'],'a_b_c'),
    (['_a', 'b.c'],'a_b.c'),
    (['a', '-b.--'],'a_b'),
    ([' ', '-b', 'c '],'b_c'),
    (['', 'b'],'b'),
])
def test_join_names(input, expected):
    assert join_names(input) == expected

#TODO
def test_get_atoms():
    pass

#TODO
def test_get_images():
    pass
