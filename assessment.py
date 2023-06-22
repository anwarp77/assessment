import sqlite3
import re
from Bio.Seq import Seq


def primers():
    """
    This funtion reads in the primers.fasta file and creates a dictionary
    that has the markers as key and the forwards and revers primers
    as value
    :return dict_primers: dictionary with marker as key and forward and reverse
    primers as value
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


def checkprimer(cp_primers, cp_read, c):
    """
    This function searches for the primer that is present in the given
    read. It first loops through the dictionary and create the variables
    for the original primers and reverse complement (RC) primers. By using
    re.search and regex it searches for the fw_primer + sequence in
    between + rv_primer in the original and RC sequence. if its a
    match it splices the sequence and the quality score (QS). Furthermore,
    the Header, sequence, QS, RC and Marker will be stored in the database
    afterwards.
    :param d_primers: primer dictionary
    :param l_reads: a read (header, sequence, +, QS)
    :param c: cursor database
    :return: None
    """
    sequence = cp_read[1]
    rev_comp_seq = Seq(sequence).reverse_complement()
    score = cp_read[3]
    RC = None

    for marker, primer in cp_primers.items():
        fw_primer = primer[0]
        rv_primer = primer[1]
        revcomp_fw_primer = Seq(fw_primer).reverse_complement()
        revcomp_rv_primer = Seq(rv_primer).reverse_complement()

        #search for the fw_primer and rv_primer in the original sequence and RC sequence
        fw_match = re.search(f"{fw_primer}.*{rv_primer}", sequence)
        rv_match = re.search(f"{fw_primer}.*{rv_primer}", str(rev_comp_seq))

        #splicing the sequence and RC sequence if it has a match
        if fw_match:
            f_index = fw_match.start() + len(fw_primer)
            r_index = fw_match.end() - len(rv_primer)
            sequence = sequence[f_index:r_index]
            score = score[f_index:r_index]
            RC = "False"
        if rv_match:
            f_index = rv_match.start() + len(revcomp_fw_primer)
            r_index = rv_match.end() - len(revcomp_rv_primer)
            sequence = str(rev_comp_seq)[f_index:r_index]
            score = score[f_index:r_index]
            RC = "True"

        # if RC is None it means the read doesn't contain a primer and will not be stored in the database
        if RC is not None:
            insert_read = """
                INSERT INTO reads VALUES (?, ?, ?, ?, ?)
            """
            c.execute(insert_read,
                      (cp_read[0], sequence, score, RC, marker))
            break


def inlezen(c):
    """
    This function reads in the reads.fastq file and check the reads one
    by one if there it contains one of the primers using the checkprimer()
    funtion.
    :param c: cursor database
    :return: None
    """
    with open("reads.fastq", "r") as reads_infile:
        lijst_read = []
        for line in reads_infile:
            line = line.strip("\n")
            lijst_read.append(line)
            # check for primers one by one
            if len(lijst_read) == 4:
                checkprimer(primers(), lijst_read, c)
                lijst_read = []


def db_queries(c):
    """
    This function prints out interesting queries
    :param c: cursor database
    :return: None
    """
    c.execute(
        "SELECT DISTINCT Marker, count(Marker) FROM reads GROUP BY Marker HAVING RC LIKE 'True' ORDER BY count(Marker) DESC")
    results = c.fetchall()
    print(
        f"Most common markers in the reverse complement sequence: \n{results}\n", )

    c.execute(
        "SELECT DISTINCT Marker, count(Marker) FROM reads GROUP BY Marker HAVING RC LIKE 'False' ORDER BY count(Marker) DESC")
    results = c.fetchall()
    print(f"Most common markers in de original sequence: \n{results}\n")

    c.execute(
        "SELECT sum(length(sequence))/(SELECT count(*) FROM reads) FROM reads")
    results = c.fetchall()
    print(f"Average sequence length: \n{results}\n")


if __name__ == '__main__':
    #connect/create the database
    conn = sqlite3.connect('markerSeq.db')
    c = conn.cursor()
    c.execute("""
        DROP TABLE IF EXISTS reads
    """)
    c.execute("""
        CREATE TABLE IF NOT EXISTS reads(
        Header VARCHAR(50) NOT NULL PRIMARY KEY,
        Sequence TEXT NOT NULL,
        QS TEXT NOT NULL,
        RC BOOLEAN NOT NULL,
        Marker VARCHAR NOT NULL
        )
    """)

    data_db = inlezen(c)
    db_queries(c)

    conn.commit()

    conn.close()
