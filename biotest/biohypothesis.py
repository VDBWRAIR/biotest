from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from hypothesis import strategies as st

##############
# Hypothesis #
##############
make_seqrec = lambda id, seq, quals: \
                SeqRecord(Seq(seq, IUPAC.ambiguous_dna), id=str(id), description='', letter_annotations={'phred_quality':quals})
rec = st.integers(min_value=1, max_value=10).flatmap(
    lambda n:
        st.builds(
            make_seqrec,
            st.integers(),
            st.text(alphabet='ATGCN', min_size=n, max_size=n),
            st.lists(st.integers(min_value=35, max_value=40), min_size=n, max_size=n)
        )
)
