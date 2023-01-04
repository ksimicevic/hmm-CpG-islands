import random

NAME = "chr19"

PATH_TO_SEQUENCE = "big-sequences/" + "chr19.fa"
PATH_TO_ISLANDS = "data/strict-examples/" + "chr19_islands.csv"

PATH_TO_END_RES = "data/even/sequences/"
PATH_TO_END_ISLANDS = "data/even/islands/"

NUM_OF_ISLANDS = 3

file = open(PATH_TO_SEQUENCE, 'r')
seq = ''.join([s.strip().capitalize() for s in file.read()])

file = open(PATH_TO_ISLANDS, 'r')
islands = [s.strip().split(',') for s in file.readlines()]

for i in range(0, len(islands), NUM_OF_ISLANDS):
    sequence = ""
    indexes = []

    for j in range(0, NUM_OF_ISLANDS):
        island = islands[i + j]

        begin, end = int(island[0]), int(island[1])
        island_len = end - begin + 1

        begin_perc = random.uniform(0, 1)
        end_perc = 1 - begin_perc

        begin_seq, end_seq = begin - int(begin_perc * island_len), end + int(end_perc * island_len)

        indexes.append((begin - begin_seq + len(sequence), end - begin_seq + len(sequence)))
        sequence += seq[begin_seq: end_seq]

    result_file = open(f"{PATH_TO_END_RES}{NAME}_{i}.txt", 'w+')
    result_file.write(sequence)

    result_island_file = open(f"{PATH_TO_END_ISLANDS}{NAME}_{i}.txt", 'w+')
    result_island_file.write('\n'.join([f"{beg},{end}" for (beg, end) in indexes]))
