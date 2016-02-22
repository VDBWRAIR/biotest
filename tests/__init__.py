try:
    import unittest2 as unittest
except:
    import unittest
import sys
if sys.version[2] != '6':
    from hypothesis import strategies as st
    st.text().example()


