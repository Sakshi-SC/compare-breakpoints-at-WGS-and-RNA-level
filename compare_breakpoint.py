####This script is designed to compare structural variant (SV) breakpoints identified at the whole genome sequencing (WGS) level using BreakDancer with fusion breakpoints identified at the RNA level using STAR-Fusion#####

def parse_position_info(file_path):
    breakpoints = []
    with open(file_path, 'r') as file:
        # Skip header line
        next(file)
        for line in file:
            # Split the line by tab or space
            chrom1_pos1, chrom2_pos2 = line.strip().split('\t')
            # Check if the line contains the header
            if chrom1_pos1.startswith('#'):
                continue  # Skip header line
            # Split chromosome and position for each breakpoint
            chrom1, pos1 = chrom1_pos1.split(':')
            chrom2, pos2 = chrom2_pos2.split(':')
            # Append to the list of breakpoints
            breakpoints.append((chrom1, int(pos1), chrom2, int(pos2)))
    return breakpoints
def compare_breakpoints(bd_breakpoints, sf_breakpoints, threshold):
    matched_pairs = []
    for bd_bp in bd_breakpoints:
        for sf_bp in sf_breakpoints:
            #if bd_bp[0] == sf_bp[0] and bd_bp[2] == sf_bp[2]:
                dist1 = abs(bd_bp[1] - sf_bp[1])
                dist2 = abs(bd_bp[3] - sf_bp[3])
                if dist1 <= threshold and dist2 <= threshold:
                    matched_pairs.append((bd_bp, sf_bp))
    return matched_pairs

# Parse breakpoint information from files
bd_breakpoints = parse_position_info('BD_output')
sf_breakpoints = parse_position_info('SF_output')

# Define threshold for position matching in base pairs
threshold = 500

# Compare breakpoints
matched_pairs = compare_breakpoints(bd_breakpoints, sf_breakpoints, threshold)

# Report results
for bd_bp, sf_bp in matched_pairs:
    print("Matched Pair:")
    print("BreakDancer Breakpoint:", bd_bp)
    print("STAR-Fusion Breakpoint:", sf_bp)
    print()


