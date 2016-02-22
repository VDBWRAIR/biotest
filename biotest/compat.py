try:
    from itertools import takewhile, imap, ifilter, izip
except ImportError:
    imap = map
    ifilter = filter
    izip = zip
