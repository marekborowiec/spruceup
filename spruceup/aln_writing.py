# Functions to write popular multiple sequence alignment file formats
# from  dict {taxon: sequence}


def print_fasta(source_dict):
    """Print fasta-formatted string from {taxon: seq} dict."""
    fasta_string = ""
    # each sequence line will have 80 characters
    n = 80
    for taxon, seq in sorted(source_dict.items()):
        # split dictionary values to a list of string, each n chars long
        seq = [seq[i : i + n] for i in range(0, len(seq), n)]
        # in case there are unwanted spaces in taxon names
        taxon = taxon.replace(" ", "_").strip("'")
        fasta_string += ">" + taxon + "\n"
        for element in seq:
            fasta_string += element + "\n"
    return fasta_string


def print_phylip(source_dict):
    """Print phylip-formatted string from {taxon: seq} dict."""
    taxa_list = list(source_dict.keys())
    no_taxa = len(taxa_list)
    # figure out the max length of a taxon for nice padding of sequences
    pad_longest_name = len(max(taxa_list, key=len)) + 3
    # get sequence length from a random value
    seq_length = len(next(iter(source_dict.values())))
    header = str(len(source_dict)) + " " + str(seq_length)
    phylip_string = header + "\n"
    for taxon, seq in sorted(source_dict.items()):
        taxon = taxon.replace(" ", "_").strip("'")
        # left-justify taxon names relative to sequences
        phylip_string += taxon.ljust(pad_longest_name, ' ') + seq + "\n"
    return phylip_string


def print_phylip_int(source_dict):
    """Print phylip interleaved-formatted string from {taxon: seq} dict."""
    taxa_list = list(source_dict.keys())
    no_taxa = len(taxa_list)
    pad_longest_name = len(max(taxa_list, key=len)) + 3
    seq_length = len(next(iter(source_dict.values())))
    header = str(len(source_dict)) + " " + str(seq_length)
    phylip_int_string = header + "\n\n"
    # this will be a list of tuples to hold taxa names and sequences
    seq_matrix = []
    # each sequence line will have 500 characters
    n = 500
    # recreate sequence matrix
    add_to_matrix = seq_matrix.append
    for taxon, seq in sorted(source_dict.items()):
        add_to_matrix((taxon, [seq[i : i + n] for i in range(0, len(seq), n)]))
    first_seq = seq_matrix[0][1]
    for index, item in enumerate(first_seq):
        for taxon, sequence in seq_matrix:
            if index == 0:
                phylip_int_string += (
                    taxon.ljust(pad_longest_name, ' ') + sequence[index] + "\n"
                )
            else:
                phylip_int_string += sequence[index] + "\n"
        phylip_int_string += "\n"
    return phylip_int_string


def print_nexus(source_dict, data_type):
    """Print nexus-formatted string from {taxon: seq} dict."""
    if data_type == "aa" or command == "translate":
        data_type = "PROTEIN"
    elif data_type == "nt":
        data_type = "DNA"
    taxa_list = list(source_dict.keys())
    no_taxa = len(taxa_list)
    pad_longest_name = len(max(taxa_list, key=len)) + 3
    seq_length = len(next(iter(source_dict.values())))
    header = str(len(source_dict)) + " " + str(seq_length)
    nexus_string = (
        "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS  NTAX="
        + str(no_taxa)
        + " NCHAR="
        + str(seq_length)
        + ";\n\tFORMAT DATATYPE="
        + data_type
        + "  GAP = - MISSING = ?;\n\tMATRIX\n"
    )
    for taxon, seq in sorted(source_dict.items()):
        taxon = taxon.replace(" ", "_").strip("'")
        nexus_string += "\t" + taxon.ljust(pad_longest_name, ' ') + seq + "\n"
    nexus_string += "\n;\n\nEND;"
    return nexus_string


def print_nexus_int(source_dict, data_type):
    """Print nexus interleaved-formatted string from {taxon: seq} dict."""
    if data_type == "aa":
        data_type = "PROTEIN"
    elif data_type == "nt":
        data_type = "DNA"
    taxa_list = list(source_dict.keys())
    no_taxa = len(taxa_list)
    pad_longest_name = len(max(taxa_list, key=len)) + 3
    seq_length = len(next(iter(source_dict.values())))
    header = str(len(source_dict)) + " " + str(seq_length)
    # this will be a list of tuples to hold taxa names and sequences
    seq_matrix = []
    nexus_int_string = (
        "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS  NTAX="
        + str(no_taxa)
        + " NCHAR="
        + str(seq_length)
        + ";\n\tFORMAT   INTERLEAVE"
        + "   DATATYPE="
        + data_type
        + "  GAP = - MISSING = ?;\n\tMATRIX\n"
    )
    n = 500
    # recreate sequence matrix
    add_to_matrix = seq_matrix.append
    for taxon, seq in sorted(source_dict.items()):
        add_to_matrix((taxon, [seq[i : i + n] for i in range(0, len(seq), n)]))
    first_seq = seq_matrix[0][1]
    for index, item in enumerate(first_seq):
        for taxon, sequence in seq_matrix:
            if index == 0:
                nexus_int_string += (
                    taxon.ljust(pad_longest_name, ' ') + sequence[index] + "\n"
                )
            else:
                nexus_int_string += sequence[index] + "\n"
        nexus_int_string += "\n"
    nexus_int_string += "\n;\n\nEND;"
    return nexus_int_string


def write_alignment_file(alignment_dict, file_format, out_file_name, data_type):
    """Write the correct format string to file."""
    with open(out_file_name, 'w') as of:
        if file_format == 'phylip':
            of.write(print_phylip(alignment_dict))
        elif file_format == 'fasta':
            of.write(print_fasta(alignment_dict))
        elif file_format == 'phylip-int':
            of.write(print_phylip_int(alignment_dict))
        elif file_format == 'nexus':
            of.write(print_nexus(alignment_dict, data_type))
        elif file_format == 'nexus-int':
            of.write(print_nexus_int(alignment_dict, data_type))

def main():
    pass

if __name__ == '__main__':
    main()
