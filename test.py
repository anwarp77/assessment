import re
from Bio.Seq import Seq
string = "AATAGAGCAGTGTCGCGCACTCTGCGTCACGATGTGTGGTCTATGACAGTTCTTTTTTTGTGTGGTTACGCGGACCACCGAGATTTCACCCGCCGTGGATATAC"
rev_com = Seq(string).reverse_complement()
f_primer = Seq("AGTGTCGCGCACT").reverse_complement()
r_primer = Seq("ACCGAGATTTCAC").reverse_complement()
result = re.search(f'{f_primer}.*{r_primer}', rev_com)
print(result)