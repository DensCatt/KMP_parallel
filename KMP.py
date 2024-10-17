import os
import argparse
import multiprocessing

parser = argparse.ArgumentParser(description="DNA Pattern Matching using KMP algorithm")
parser.add_argument('--DNA', type=str, default=None, help='DNA sequence in FASTA format', required=True)
parser.add_argument("--patterns", type=str, default=None,
                    help="Patterns to match on the DNA sequence (FASTA format)", required=True)
parser.add_argument('--output', type=str, help="txt output file (default: patterns_matches.txt)",
                    default="patterns_matches.txt")
parser.add_argument('--num_cores', type=str, default=4, help='# cores for parallel execution (default: 4)')

args2 = parser.parse_args()


# Exception to handle empty sequences in the patterns file
class EmptySequenceException(Exception):
    def __init__(self, msg):
        super().__init__(msg)
        self.msg = msg


# Function to index the pattern
def create_table(pattern):
    j = 0
    i = 1
    table = [0] * len(pattern)  # initialize the table with zeros (it is actually a list of numbers)

    while i < len(pattern):
        # in case of match with the characters change the correspondent value and move forward
        if pattern[i] == pattern[j]:
            j += 1
            table[i] = j
            i += 1
        # when there is not a match check the value of j
        else:
            if j != 0:
                j = table[j - 1]
            else:
                table[i] = 0  # no matches
                i += 1
    return table


# Function to perform the Knuth-Morris-Pratt search of a pattern on a sequence
def KMPsearch(pattern, string):
    pat_len = len(pattern)
    # If the sequence is empty raises an exception
    if pat_len == 0:
        raise EmptySequenceException("Empty pattern")

    str_len = len(string)
    table = create_table(pattern)  # Indexing the pattern
    indexes = []

    i = 0
    j = 0
    while i < str_len:
        if pattern[j] == string[i]:
            j += 1
            i += 1
            if j == pat_len:  # the end of the pattern
                indexes.append(i - j)  # add the position of the match to the list
                j = table[j - 1]

        else:
            if j != 0:
                j = table[j - 1]
            else:
                i += 1

    return indexes


# Function to read a FASTA file
def read_fasta(filename):
    # The check on the file has been already done
    # read the file and save the sequence in a string
    print("Reading the DNA sequence file...")
    with open(filename) as f:
        f.readline()
        DNA_seq = "".join([line.strip() for line in f])
        
        if not DNA_seq:
            raise EmptySequenceException("Empty DNA sequence file")
        if ">" in DNA_seq:
            raise ValueError("--DNA file must contain just one sequence")
    return DNA_seq


# Function to read a FASTA file with multiple sequences
def read_multifasta(filename):
    # The check on the file has been already done
    # read the file and save the sequences in a list of strings
    print("Reading the pattern file...")
    with open(filename) as f:
        f.readline()
        pattern = ""
        pattern_list = []
        # NOTE: a more efficient method would be:
        # pattern_list = [line.strip() for line in f if not line.startswith(">")]
        # but this does not work properly in case of long multi-line patterns...
        for line in f:
            if not line.startswith(">"):
                pattern += line.strip()
            else:
                pattern_list.append(pattern)
                pattern = ""
        pattern_list.append(pattern)
    return pattern_list


# Function to be executed by each process
def worker(args):
    try:
        pattern, DNA_sequence, lock, output_file, empty_pat = args
        matches_idx = KMPsearch(pattern, DNA_sequence)
        with lock:  # manages the access to the file

            # write on the output file the position of each match
            with open(output_file, 'a') as f:
                if len(matches_idx) == 0:
                    f.write(f"Pattern: {pattern}, no matches found\n")
                for idx in matches_idx:
                    f.write(f"Pattern: {pattern}, Match at index: {idx}\n")

    except EmptySequenceException as e:
        empty_pat.value += 1
        print(f"WARNING: {e.msg}")


def main():
    errors = []
    # Handling the possible input errors BEFORE the code execution
    # Saving them on a list to show them all at once
    if not args2.output.endswith(".txt"):
        errors.append("ERROR: Output file must be a text file")

    if not os.path.exists(args2.DNA):
        errors.append(f"ERROR: DNA file '{args2.DNA}' not found")
    else:
        with open(args2.DNA) as f:
            if not f.readline().startswith(">"):
                errors.append(f"ERROR: DNA file '{args2.DNA}' must be in FASTA format")

    if not os.path.exists(args2.patterns):
        errors.append(f"ERROR: Pattern file '{args2.patterns}' not found")
    else:
        with open(args2.patterns) as f:
            if not f.readline().startswith(">"):
                errors.append(f"ERROR: DNA file '{args2.patterns}' must be in FASTA format")

    if errors:
        for error in errors:
            print(error)
        print()
        return  # quits

    DNA_sequence = read_fasta(args2.DNA)
    patterns = read_multifasta(args2.patterns)
    cores = int(args2.num_cores)
    output_file = args2.output

    if cores > multiprocessing.cpu_count():
        raise ValueError(f"Available cores are less than {cores}")

    # manager and lock to safely access the output file
    manager = multiprocessing.Manager()
    lock = manager.Lock()

    empty_pat = manager.Value('i', 0)

    with open(output_file, 'w') as f:
        f.write(">TEXT FILE CONTAINING THE MATCHES FOUND FOR EACH PATTERN:\n\n")  # header

    # Create a pool of processes to perform the algorithm in parallel
    print(f"Using {cores} cores for processing...")
    pool = multiprocessing.Pool(processes=cores)
    args = [(pattern, DNA_sequence, lock, output_file, empty_pat) for pattern in patterns]
    pool.map(worker, args)
    pool.close()
    pool.join()
    print("Done! Open the output file to see the results\n")

    if empty_pat.value > 0:
        if empty_pat.value == 1:
            print(f"WARNING: There was 1 empty sequence in the patterns file\n")
        else:
            print(f"WARNING: There were {empty_pat.value} empty sequences in the patterns file\n")


if __name__ == "__main__":
    try:
        main()
    except EmptySequenceException as e:
        print(f"ERROR: {e.msg}")
    except ValueError as e:
        print(f"ERROR: {e}")
