import os.path, sys.argv

vcf_filename = sys.argv[1]
prefix = sys.argv[2]
suffix = sys.argv[3]

results = []
with open(vcf_filename, 'r') as vcf_file:
    for line in vcf_file:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        chrom, start, info = fields[0], int(fields[1]), fields[8].split(':')
        end = int(start) + 1
        samples = [sample.split(':') for sample in fields[9:]]
        values = []
        for sample in samples:
            d  = {}
            for i in range(len(sample)):
                d[info[i]] = sample[i]
            if 'AD' in d:
                value = int(d['AD'])
            else:
                value = 0
            values.append(value)
        results.append([chrom, start, end] + values)
print(len(results), 'lines processed;')
nsamples = len(results[0][3:])
for n in range(nsamples):
    sample_name = 's' + '{:02d}'.format(n+1) + '_' + suffix
    output_filename = os.path.join(os.path.dirname(vcf_filename), prefix, '_' + sample_name + '_' + suffix + '.bedgraph')
    with open(output_filename, 'w') as output_file:
        count = 0
        for result in results:
            value = result[n+3]
            if value > 0:
                count += 1
                out = result[:3] + [value]
                output_file.write('\t'.join([str(f) for f in out]) + '\n')
        print('Sample', (n+1), '=', count, 'non-zero values;')
print('Done!')