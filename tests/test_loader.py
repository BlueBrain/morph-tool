import mock
from mock import patch

import morph_tool.loader as tested


def test_ensure_startswith_point():
    assert (
        tested._ensure_startswith_point(".ext") ==
        ".ext")
    assert (
        tested._ensure_startswith_point("ext") ==
        ".ext")


@patch('morphio.Morphology')
def test_loader(f_mock):
    f_mock.configure_mock(side_effect=lambda *args: object())
    loader = tested.MorphLoader('/dir', file_ext='abc', cache_size=1)
    morph1 = loader.get('test')
    # should get cached object now
    assert (
        loader.get('test') is
        morph1)
    # options are different => should not get cached object
    assert (
        loader.get('test', options=42) is not
        morph1)
    # first cached object was evicted from the cache
    assert (
        loader.get('test') is not
        morph1)
    f_mock.assert_has_calls([
        mock.call('/dir/test.abc'),
        mock.call('/dir/test.abc', 42),
        mock.call('/dir/test.abc'),
    ])

@patch('morphio.Morphology')
def test_loader_no_cache(f_mock):
    f_mock.configure_mock(side_effect=lambda *args: object())
    loader = tested.MorphLoader('/dir', file_ext='abc', cache_size=0)
    loader.get('test')
    loader.get('test')
    assert f_mock.call_count == 2
