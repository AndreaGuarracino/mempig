import sys

def load_info(info_file):
    info = {}
    with open(info_file, 'r') as f:
        for line in f:
            name, length = line.strip().split()
            info[name] = int(length)
    return info

def convert_to_paf(input_file, output_file, target_info):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split()
            query_name = fields[0]
            query_start = int(fields[1])  # 0-based
            query_end = int(fields[2])
            
            mem_length = query_end - query_start
            
            # Process all target fields (column 6 and beyond)
            target_fields = fields[5:]
            for target_full in target_fields:
                target_parts = target_full.split(':')
                target_name = ':'.join(target_parts[:-2])
                strand = target_parts[-2]
                target_start = int(target_parts[-1])
                target_end = target_start + mem_length

                query_length = 150  # Fixed query length
                target_length = target_info.get(target_name, 0)
                
                # PAF fields
                paf_fields = [
                    query_name,
                    str(query_length),
                    str(query_start),
                    str(query_end),
                    strand,
                    target_name,
                    str(target_length),
                    str(target_start),
                    str(target_end),
                    str(mem_length),        # Number of matching bases
                    str(mem_length),        # Alignment block length
                    '255',                  # Mapping quality (using 255 as a placeholder)
                    f'cg:Z:{mem_length}='   # CIGAR string
                ]
                
                outfile.write('\t'.join(paf_fields) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <target_info_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    target_info_file = sys.argv[2]
    output_file = sys.argv[3]

    target_info = load_info(target_info_file)
    convert_to_paf(input_file, output_file, target_info)
