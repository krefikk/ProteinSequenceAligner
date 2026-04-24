from Bio.Align import substitution_matrices

def slice_sequence(seq, start, end):
    idx_start = max(0, start - 1) if start else 0
    idx_end = end if end else len(seq)
    return seq[idx_start:idx_end]

def get_blosum_score(matrix, aa1, aa2):
    try:
        return matrix[aa1, aa2]
    except KeyError:
        return matrix[aa2, aa1]

def align_dp(seq1, seq2, mode='global', gap_penalty=10.0, seq1_bounds=(None, None), seq2_bounds=(None, None)):
    s1 = slice_sequence(seq1, seq1_bounds[0], seq1_bounds[1])
    s2 = slice_sequence(seq2, seq2_bounds[0], seq2_bounds[1])
    
    blosum62 = substitution_matrices.load("BLOSUM62")
    gap = abs(gap_penalty)
    
    n = len(s1)
    m = len(s2)
    
    score_matrix = [[0 for j in range(m + 1)] for i in range(n + 1)]
    traceback = [[None for j in range(m + 1)] for i in range(n + 1)]
    
    if mode == 'global':
        for i in range(1, n + 1):
            score_matrix[i][0] = -i * gap
            traceback[i][0] = 'up'
        for j in range(1, m + 1):
            score_matrix[0][j] = -j * gap
            traceback[0][j] = 'left'
            
    max_local_score = 0
    max_local_pos = (0, 0)
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = get_blosum_score(blosum62, s1[i-1], s2[j-1])
            
            diag = score_matrix[i-1][j-1] + match_score
            up   = score_matrix[i-1][j] - gap
            left = score_matrix[i][j-1] - gap
            
            if mode == 'global':
                scores = [diag, up, left]
                best = max(scores)
                score_matrix[i][j] = best
                
                if best == diag: traceback[i][j] = 'diag'
                elif best == up: traceback[i][j] = 'up'
                else:            traceback[i][j] = 'left'
                
            elif mode == 'local':
                scores = [0, diag, up, left]
                best = max(scores)
                score_matrix[i][j] = best
                
                if best == 0:    traceback[i][j] = 'stop'
                elif best == diag: traceback[i][j] = 'diag'
                elif best == up: traceback[i][j] = 'up'
                else:            traceback[i][j] = 'left'
                
                if best > max_local_score:
                    max_local_score = best
                    max_local_pos = (i, j)

    align1, align2, lines = "", "", ""
    
    if mode == 'global':
        i, j = n, m
        final_score = score_matrix[n][m]
    else:
        i, j = max_local_pos
        final_score = max_local_score

    path = [(i, j)]

    while True:
        if mode == 'global' and i == 0 and j == 0:
            break
        if mode == 'local' and (traceback[i][j] == 'stop' or score_matrix[i][j] == 0):
            break
            
        current_dir = traceback[i][j]
        
        if current_dir == 'diag':
            align1 += s1[i-1]
            align2 += s2[j-1]
            lines += "|" if s1[i-1] == s2[j-1] else "."
            i -= 1
            j -= 1
        elif current_dir == 'up':
            align1 += s1[i-1]
            align2 += "-"
            lines += " "
            i -= 1
        elif current_dir == 'left':
            align1 += "-"
            align2 += s2[j-1]
            lines += " "
            j -= 1
            
        path.append((i, j))

    align1 = align1[::-1]
    align2 = align2[::-1]
    lines = lines[::-1]

    return final_score, align1, lines, align2, score_matrix, path

def format_alignment_text(seq1, lines, seq2, chunk_size=80):
    result = []
    for i in range(0, len(seq1), chunk_size):
        result.append(seq1[i:i+chunk_size])
        result.append(lines[i:i+chunk_size])
        result.append(seq2[i:i+chunk_size])
        result.append("-" * chunk_size)
    return "\n".join(result)

def parse_fasta(filepath):
    sequence = ""
    try:
        with open(filepath, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                sequence += line
        return sequence.upper()
    except Exception as e:
        return None

# For testing without GUI
if __name__ == "__main__":
    
    # MSH2 Homo Sapiens
    seq_a = "MAVQPKETLQLESAAEVGFVRFFQGMPEKPTTTVRLFDRGDFYTAHGEDALLAAREVFKTQGVIKYMGPAGAKNLQSVVLSKMNFESFVKDLLLVRQYRVEVYKNRAGNKASKENDWYLAYKASPGNLSQFEDILFGNNDMSASIGVVGVKMSAVDGQRQVGVGYVDSIQRKLGLCEFPDNDQFSNLEALLIQIGPKECVLPGGETAGDMGKLRQIIQRGGILITERKKADFSTKDIYQDLNRLLKGKKGEQMNSAVLPEMENQVAVSSLSAVIKFLELLSDDSNFGQFELTTFDFSQYMKLDIAAVRALNLFQGSVEDTTGSQSLAALLNKCKTPQGQRLVNQWIKQPLMDKNRIEERLNLVEAFVEDAELRQTLQEDLLRRFPDLNRLAKKFQRQAANLQDCYRLYQGINQLPNVIQALEKHEGKHQKLLLAVFVTPLTDLRSDFSKFQEMIETTLDMDQVENHEFLVKPSFDPNLSELREIMNDLEKKMQSTLISAARDLGLDPGKQIKLDSSAQFGYYFRVTCKEEKVLRNNKNFSTVDIQKNGVKFTNSKLTSLNEEYTKNKTEYEEAQDAIVKEIVNISSGYVEPMQTLNDVLAQLDAVVSFAHVSNGAPVPYVRPAILEKGQGRIILKASRHACVEVQDEIAFIPNDVYFEKDKQMFHIITGPNMGGKSTYIRQTGVIVLMAQIGCFVPCESAEVSIVDCILARVGAGDSQLKGVSTFMAEMLETASILRSATKDSLIIIDELGRGTSTYDGFGLAWAISEYIATKIGAFCMFATHFHELTALANQIPTVNNLHVTALTTEETLTMLYQVKKGVCDQSFGIHVAELANFPKHVIECAKQKALELEEFQYIGESQGYDIMEPAAKKCYLEREQGEKIIQEFLSKVKQMPFTEMSEENITIKLKQLKAEVIAKNNSFVNEIISRIKVTT"
    
    # MLH1 Homo Sapiens
    seq_b = "MSFVAGVIRRLDETVVNRIAAGEVIQRPANAIKEMIENCLDAKSTSIQVIVKEGGLKLIQIQDNGTGIRKEDLDIVCERFTTSKLQSFEDLASISTYGFRGEALASISHVAHVTITTKTADGKCAYRASYSDGKLKAPPKPCAGNQGTQITVEDLFYNIATRRKALKNPSEEYGKILEVVGRYSVHNAGISFSVKKQGETVADVRTLPNASTVDNIRSIFGNAVSRELIEIGCEDKTLAFKMNGYISNANYSVKKCIFLLFINHRLVESTSLRKAIETVYAAYLPKNTHPFLYLSLEISPQNVDVNVHPTKHEVHFLHEESILERVQQHIESKLLGSNSSRMYFTQTLLPGLAGPSGEMVKSTTSLTSSSTSGSSDKVYAHQMVRTDSREQKLDAFLQPLSKPLSSQPQAIVTEDKTDISSGRARQQDEEMLELPAPAEVAAKNQSLEGDTTKGTSEMSEKRGPTSSNPRKRHREDSDVEMVEDDSRKEMTAACTPRRRIINLTSVLSLQEEINEQGHEVLREMLHNHSFVGCVNPQWALAQHQTKLYLLNTTKLSEELFYQILIYDFANFGVLRLSEPAPLFDLAMLALDSPESGWTEEDGPKEGLAEYIVEFLKKKAEMLADYFSLEIDEEGNLIGLPLLIDNYVPPLEGLPIFILRLATEVNWDEEKECFESLSKECAMFYSIRKQYISEESTLSGQQSEVPGSIPNSWKWTVEHIVYKALRSHILPPKHFTEDGNILQLANLPDLYKVFERC"
    
    read_fasta = int(input("\nEnter 1 to read sequences from FASTA files, 0 to use hardcoded sequences: "))
    if read_fasta == 1:
        seq_a = parse_fasta(input("Enter path to Sequence A FASTA file: "))
        seq_b = parse_fasta(input("Enter path to Sequence B FASTA file: "))
    else:
        print("Using hardcoded sequences.")
        
    set_bounds_seq1 = input("\nSet bounds for Sequence A (e.g., 5-30) or press Enter to use full sequence: ")
    set_bounds_seq2 = input("Set bounds for Sequence B (e.g., 10-50) or press Enter to use full sequence: ")
    start_a, end_a = (None, None)
    start_b, end_b = (None, None)
    
    if set_bounds_seq1:
        try:
            start_a, end_a = map(int, set_bounds_seq1.split('-'))
        except ValueError:
            print("Invalid bounds format for Sequence A. Using full sequence.")
    if set_bounds_seq2:
        try:
            start_b, end_b = map(int, set_bounds_seq2.split('-'))
        except ValueError:
            print("Invalid bounds format for Sequence B. Using full sequence.")
    
    try:
        gap_penalty = float(input("\nEnter gap penalty: "))
    except ValueError:
        print("Invalid gap penalty. Using default value of 10.")
        gap_penalty = 10.0
    
    # Global alignment test
    print("\nGlobal Alignment:")
    score, align1, lines, align2, matrix, path = align_dp(seq_a, seq_b, mode='global', gap_penalty=gap_penalty, seq1_bounds=(start_a, end_a), seq2_bounds=(start_b, end_b))
    print(f"Score: {score}")
    print(format_alignment_text(align1, lines, align2))
    
    # Local alignment test
    print("\nLocal Alignment:")
    score, align1, lines, align2, matrix, path = align_dp(seq_a, seq_b, mode='local', gap_penalty=gap_penalty, seq1_bounds=(start_a, end_a), seq2_bounds=(start_b, end_b))
    print(f"Score: {score}")
    print(format_alignment_text(align1, lines, align2))