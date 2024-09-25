import io

from gpsea.util import open_text_io_handle_for_reading


def test_bytes():
    buf = io.BytesIO(b"some binary data: 123")
    with open_text_io_handle_for_reading(buf) as fh:
        data = fh.read()

    assert data == 'some binary data: 123'


def test_text():
    buf = io.StringIO("some text data: 123")
    with open_text_io_handle_for_reading(buf) as fh:
        data = fh.read()

    assert data == 'some text data: 123'
