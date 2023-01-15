import os
import errno
import argparse

from Bio import SeqIO
from io import StringIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation

parser = argparse.ArgumentParser(description='Visualise a sequence')
parser.add_argument('sequence_path', metavar='sequence_path', action='store', type=str, help='Path to the sequence')
parser.add_argument('real_islands_path', metavar='real_islands_path', action='store', type=str, help='Path to actual island indexes')
parser.add_argument('predicted_islands_path', metavar='predicted_islands_path', action='store', type=str, help='Path to predicted island indexes')
parser.add_argument('--output', metavar='output', action='store', type=str, help='Output file')

args = parser.parse_args()

output = 'comparison'

if args.output:
    output = args.output

raw_sequence_path = fr'{args.sequence_path}'
real_islands_path = fr'{args.real_islands_path}'
predicted_islands_path = fr'{args.predicted_islands_path}'

if not os.path.exists(raw_sequence_path):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), raw_sequence_path)

if not os.path.exists(real_islands_path):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), real_islands_path)

if not os.path.exists(predicted_islands_path):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), predicted_islands_path)


with open(raw_sequence_path, 'r') as file:
    raw_sequence = file.read()
    raw_sequence = ">seq\n" + raw_sequence

record = SeqIO.read(StringIO(raw_sequence), 'fasta')

existing_subsequences = []
rest_subsequences = []
start = 0
with open(real_islands_path, 'r') as file:
    for line in file:
        s_, e_ = line.strip().split(',')
        s_, e_ = int(s_), int(e_)
        existing_subsequences.append((s_, e_))
        if s_ > start:
            rest_subsequences.append((start, s_ - 1))
        start = e_ + 1

        feature = SeqFeature(FeatureLocation(s_, e_), type='Actual')
        record.features.append(feature)


if start < len(record):
    rest_subsequences.append((start, len(record) - 1))

for s, e in rest_subsequences:
    feature = SeqFeature(FeatureLocation(s, e), type='Rest')
    record.features.append(feature)


gdd = GenomeDiagram.Diagram("CpG Islands")
gd_track = gdd.new_track(2, name="Actual islands", greytrack=False)
gd_fs = gd_track.new_set()
for feature in record.features:
    if feature.type == "Actual":
        gd_fs.add_feature(feature, sigil="BOX", color="red", label=True, label_size=12)
    else:
        gd_fs.add_feature(feature, sigil="BOX", color="green", label=False, label_size=12)


record.features = []
gd_track2 = gdd.new_track(1, name="Predicted islands", greytrack=False)
gd_fs2 = gd_track2.new_set()

existing_subsequences = []
rest_subsequences = []
start = 0
with open(predicted_islands_path, 'r') as file:
    for line in file:
        s_, e_ = line.strip().split(',')
        s_, e_ = int(s_), int(e_)
        existing_subsequences.append((s_, e_))
        if s_ > start:
            rest_subsequences.append((start, s_ - 1))
        start = e_ + 1

        feature = SeqFeature(FeatureLocation(s_, e_), type='Predicted')
        record.features.append(feature)


if start < len(record):
    rest_subsequences.append((start, len(record) - 1))

for s, e in rest_subsequences:
    feature = SeqFeature(FeatureLocation(s, e), type='Rest')
    record.features.append(feature)
    

for feature in record.features:
    if feature.type == "Predicted":
        gd_fs2.add_feature(feature, sigil="BOX", color="blue", label=True, label_size=12)
    else:
        gd_fs2.add_feature(feature, sigil="BOX", color="green", label=False, label_size=12)

gdd.draw(format="linear", pagesize='A4', fragments=1, start=0, end=len(record))

# Check if directory images exists
if not os.path.exists('images'):
    os.makedirs('images')

gdd.write(f'images/{output}.png', 'PNG')

