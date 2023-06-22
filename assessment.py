import sqlite3
import re
from Bio.Seq import Seq

def primers():
    """
    
    """
    with open("primers.fasta", "r") as infile:
        dict_primers = {}
        marker = ""
        for line in infile:
            line = line.strip("\n")
            if line.startswith(">"):
                marker = line.strip("^>*_[RF]")
                dict_primers.setdefault(marker, [])
            else:
                dict_primers[marker].append(line)
    return dict_primers


def checkprimer(d_primers, l_reads, c):
    sequence = l_reads[1]
    rev_comp_seq = Seq(sequence).reverse_complement()
    score = l_reads[3]
    RC = None

    for marker, primer in d_primers.items():
        forward_primer = primer[0]
        reverse_primer = primer[1]

        forward_match = re.search(f"{forward_primer}.*{reverse_primer}", sequence)
        reverse_match = re.search(f"{forward_primer}.*{reverse_primer}", str(rev_comp_seq))

        if forward_match:
            f_index = forward_match.start() + len(forward_primer)
            r_index = forward_match.end() - len(reverse_primer)
            sequence = sequence[f_index:r_index]
            score = score[f_index:r_index]
            RC = "False"
        if reverse_match:
            f_index = reverse_match.start() + len(forward_primer)
            r_index = reverse_match.end() - len(reverse_primer)
            sequence = sequence[f_index:r_index]
            score = score[f_index:r_index]
            RC = "True"

        if RC is not None:
            insert_read = """
                INSERT INTO reads VALUES (?, ?, ?, ?, ?)
            """
            c.execute(insert_read, (l_reads[0], sequence, score, RC, marker))
            break
    
def inlezen(c):
    with open("reads.fastq", "r") as reads_infile:
        lijst_read = []
        data_db = []
        for line in reads_infile:
            line = line.strip("\n")
            lijst_read.append(line)
            if len(lijst_read) == 4:
                checkprimer(primers(), lijst_read, c)
                lijst_read = []
    return data_db


def db_queries(c):
    c.execute("SELECT DISTINCT Marker, count(Marker) FROM reads GROUP BY Marker HAVING RC LIKE 'True' ORDER BY count(Marker) DESC")
    results = c.fetchall()
    print(f"Most common markers in the reverse complement sequence: \n{results}\n", )

    c.execute("SELECT DISTINCT Marker, count(Marker) FROM reads GROUP BY Marker HAVING RC LIKE 'False' ORDER BY count(Marker) DESC")
    results = c.fetchall()
    print(f"Most common markers in de original sequence: \n{results}\n")

    c.execute("SELECT sum(length(sequence))/(SELECT count(*) FROM reads) FROM reads")
    results = c.fetchall()
    print(f"Average sequence length: \n{results}\n")




if __name__ == '__main__':
    conn = sqlite3.connect('markerSeq.db')
    c = conn.cursor()
    c.execute("""
        DROP TABLE reads
    """)
    c.execute("""
        CREATE TABLE IF NOT EXISTS reads(
        Header VARCHAR(50) NOT NULL PRIMARY KEY,
        Sequence TEXT NULL,
        QS TEXT NULL,
        RC BOOLEAN NULL,
        Marker VARCHAR
        )
    """)
    
    data_db = inlezen(c)
    db_queries(c)

    conn.commit()

    conn.close()
