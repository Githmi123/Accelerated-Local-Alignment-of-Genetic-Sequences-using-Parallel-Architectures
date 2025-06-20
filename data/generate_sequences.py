import random

def generate_sequence(length):
    return ''.join(random.choices('ACGT', k=length)) # random.choices('ACGT', k=length) will output a list like ['A', 'G', 'C', 'T', 'T', 'A', 'C', 'G', 'G', 'T']

with open("DNASequences.txt", "w") as file:
    for i in range(10000):
        seq1 = generate_sequence(200)
        seq2 = generate_sequence(200)
        file.write(f"{seq1},{seq2}\n")