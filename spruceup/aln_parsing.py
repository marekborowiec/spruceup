# Functions to parse popular multiple sequence alignment file formats
# into dict {taxon: sequence}.


import re


def read_in_file(in_file):
    """Read input file."""
    with open(in_file, 'r') as f:
        file_lines = f.read().rstrip('\r\n')
    return file_lines


def fasta_parse(in_file_lines):
    """Parse names and sequences in sequential fasta files using regex."""
    matches = re.finditer(r"^>(.*[^$])([^>]*)", in_file_lines, re.MULTILINE)
    records = {}
    for match in matches:
        name_match = match.group(1).replace("\n", "")
        seq_match = match.group(2).replace("\n", "").upper()
        seq_match = translate_ambiguous(seq_match)
        records[name_match] = seq_match
    return records


def phylip_parse(in_file_lines):
    """Parse names and sequences in sequential phylip files using regex."""
    matches = re.finditer(
        r"^(\s+)?(\S+)\s+([A-Za-z*?.{}-]+)", in_file_lines, re.MULTILINE
    )
    records = {}
    for match in matches:
        name_match = match.group(2).replace("\n", "")
        seq_match = match.group(3).replace("\n", "").upper()
        seq_match = translate_ambiguous(seq_match)
        records[name_match] = seq_match
    return records


def phylip_interleaved_parse(in_file_lines):
    """Parse names and sequences in interleaved phylip files using regex."""
    tax_chars_matches = re.finditer(
        r"^(\s+)?([0-9]+)[ \t]+([0-9]+)", in_file_lines, re.MULTILINE
    )
    name_matches = re.finditer(
        r"^(\s+)?(\S+)[ \t]+[A-Za-z*?.{}-]+", in_file_lines, re.MULTILINE
    )
    seq_matches = re.finditer(
        r"(^(\s+)?\S+[ \t]+|^)([A-Za-z*?.{}-]+)$", in_file_lines, re.MULTILINE
    )
    # get number of taxa and chars
    for match in tax_chars_matches:
        tax_match = match.group(2)
        chars_match = match.group(3)
    # initiate lists for taxa names and sequence strings on separate lines
    taxa = []
    sequences = []
    # initiate a dictionary for the name:sequence records
    records = {}
    # initiate a counter to keep track of sequences strung together
    # from separate lines
    counter = 0
    for match in name_matches:
        name_match = match.group(2).replace("\n", "")
        taxa.append(name_match)
    for match in seq_matches:
        seq_match = match.group(3).replace("\n", "").upper()
        seq_match = translate_ambiguous(seq_match)
        sequences.append(seq_match)
    # try parsing PHYLUCE-style interleaved phylip
    if len(taxa) != int(tax_match):
        taxa = []
        sequences = []
        matches = re.finditer(
            r"(^(\s+)?(\S+)( ){2,}|^\s+)([ A-Za-z*?.{}-]+)",
            in_file_lines,
            re.MULTILINE,
        )
        for match in matches:
            try:
                name_match = match.group(3).replace("\n", "")
                taxa.append(name_match)
            except AttributeError:
                pass
            seq_match = match.group(5).replace("\n", "").upper()
            seq_match = "".join(seq_match.split())
            seq_match = translate_ambiguous(seq_match)
            sequences.append(seq_match)
    for taxon_no in range(len(taxa)):
        sequence = ""
        for index in range(counter, len(sequences), len(taxa)):
            sequence += sequences[index]
        records[taxa[taxon_no]] = sequence
        counter += 1
    return records


def nexus_parse(in_file_lines):
    """Parse names and sequences in sequential nexus files using regex."""
    # find the matrix block
    matches = re.finditer(
        r"(\s+)?(MATRIX\n|matrix\n|MATRIX\r\n|matrix\r\n)(.*?;)",
        in_file_lines,
        re.DOTALL,
    )
    records = {}
    # get names and sequences from the matrix block
    for match in matches:
        matrix_match = match.group(3)
        seq_matches = re.finditer(
            r"^(\s+)?[']?(\S+\s\S+|\S+)[']?\s+([A-Za-z*?.{}-]+)($|\s+\[[0-9]+\]$)",
            matrix_match,
            re.MULTILINE,
        )
        for match in seq_matches:
            name_match = match.group(2).replace("\n", "")
            seq_match = match.group(3).replace("\n", "").upper()
            seq_match = translate_ambiguous(seq_match)
            records[name_match] = seq_match
    return records


def nexus_interleaved_parse(in_file_lines):
    """Parse names and sequences in sequential nexus files using regex."""
    # find the matrix block
    matches = re.finditer(
        r"(\s+)?(MATRIX\n|matrix\n|MATRIX\r\n|matrix\r\n)(.*?;)",
        in_file_lines,
        re.DOTALL,
    )
    # initiate lists for taxa names and sequence strings on separate lines
    taxa = []
    sequences = []
    # initiate a dictionary for the name:sequence records
    records = {}
    for match in matches:
        matrix_match = match.group(3)
        # get names and sequences from the matrix block
        seq_matches = re.finditer(
            r"^(\s+)?[']?(\S+\s\S+|\S+)[']?\s+([A-Za-z*?.{}-]+)($|\s+\[[0-9]+\]$)",
            matrix_match,
            re.MULTILINE,
        )
        for match in seq_matches:
            name_match = match.group(2)
            if name_match not in taxa:
                taxa.append(name_match)
            seq_match = match.group(3)
            sequences.append(seq_match)
    # initiate a counter to keep track of sequences strung together
    # from separate lines
    counter = 0
    for taxon_no in range(len(taxa)):
        full_length_sequence = "".join(
            [
                sequences[index]
                for index in range(counter, len(sequences), len(taxa))
            ]
        )
        records[taxa[taxon_no]] = (
            translate_ambiguous(full_length_sequence).replace("\n", "").upper()
        )
        counter += 1
    return records


def translate_ambiguous(seq):
    """Translate ambiguous characters from curly bracket format
    to single letter format and remove spaces from sequences.
    """
    seq = seq.replace('{GT}', 'K')
    seq = seq.replace('{AC}', 'M')
    seq = seq.replace('{AG}', 'R')
    seq = seq.replace('{CT}', 'Y')
    seq = seq.replace('{CG}', 'S')
    seq = seq.replace('{AT}', 'W')
    seq = seq.replace('{CGT}', 'B')
    seq = seq.replace('{ACG}', 'V')
    seq = seq.replace('{ACT}', 'H')
    seq = seq.replace('{AGT}', 'D')
    seq = seq.replace('{GATC}', 'N')
    seq = seq.replace(' ', '')
    return seq


def parse_alignment(alignment_file_name, in_format):
    """Get alignments parsed as a list of tuples (alignment name, {taxon : sequence})."""
    alignment_file_lines = read_in_file(alignment_file_name)
    if in_format == 'fasta':
        parsed_aln = fasta_parse(alignment_file_lines)
    elif in_format == 'phylip':
        parsed_aln = phylip_parse(alignment_file_lines)
    elif in_format == 'phylip-int':
        parsed_aln = phylip_interleaved_parse(alignment_file_lines)
    elif in_format == 'nexus':
        parsed_aln = nexus_parse(alignment_file_lines)
    elif in_format == 'nexus-int':
        parsed_aln = nexus_interleaved_parse(alignment_file_lines)
    return (alignment_file_name, parsed_aln)

def main():
    pass

if __name__ == '__main__':
    main()