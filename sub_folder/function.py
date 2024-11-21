

def trim_sequence(sequence, trim_5prime=True, trim_3prime=True, ran5_prime=100, read_len=800):
    """シーケンスの5'と3'をトリミング"""
    if trim_5prime:
        n_pos_5prime = sequence[0:ran5_prime].rfind("N")
        start_pos = n_pos_5prime + 10 if n_pos_5prime != -1 else 0
    else:
        start_pos = 0

    if trim_3prime:
        n_pos_3prime = sequence[read_len:].find("N")
        end_pos = n_pos_3prime + read_len - 10 if n_pos_3prime != -1 else len(sequence)
    else:
        end_pos = len(sequence)

    return sequence[start_pos:end_pos]

def align_reads_to_reference(reads, reference_seq):
    alignments = []
    for read in reads:
        alignments.append(pairwise2.align.globalms(read, reference_seq, 5, -4, -10, -0.5, one_alignment_only=True)[0])
    return alignments

def print_alignment_with_indices(alignment, line_length=60):
    """アラインメント結果をインデックスとともに出力し、マッチしない部分を強調表示"""
    seq1_aligned, seq2_aligned, score, begin, end = alignment
    alignment_length = len(seq1_aligned)
    seq1_pos = 1
    seq2_pos = 1

    result = ""
    for i in range(0, alignment_length, line_length):
        seq1_line = seq1_aligned[i:i+line_length]
        seq2_line = seq2_aligned[i:i+line_length]
        matches = ''.join(['|' if base1 == base2 else ' ' for base1, base2 in zip(seq1_line, seq2_line)])

        # インデックスを計算
        seq1_start_index = seq1_pos
        seq2_start_index = seq2_pos
        seq1_end_index = seq1_start_index + len(seq1_line) - 1
        seq2_end_index = seq2_start_index + len(seq2_line) - 1

        # 色付けされたシーケンスを作成
        seq1_colored = ''.join([f"<span style='color:red'>{base}</span>" if base1 != base2 else base for base1, base2, base in zip(seq1_line, seq2_line, seq1_line)])
        seq2_colored = ''.join([f"<span style='color:red'>{base}</span>" if base1 != base2 else base for base1, base2, base in zip(seq1_line, seq2_line, seq2_line)])

        # 結果を追加
        result += f"<pre>read_{seq1_start_index:>6} > {seq1_colored:<{line_length}} < {seq1_end_index:>6}</pre>"
        #result += f"<pre>{'':>6}   {matches:<{line_length}}   {'':>6}</pre>"
        result += f"<pre>target{seq2_start_index:>6} > {seq2_colored:<{line_length}} < {seq2_end_index:>6}</pre>"
        result += "\n"

        seq1_pos = seq1_end_index + 1
        seq2_pos = seq2_end_index + 1

    return result





def format_alignment(seq1, seq2, read_name="read", target_name="target", line_length=60):
    """アラインメント結果をインデックスとともに出力し、配列を揃えて表示"""
    formatted_result = ""

    seq1_pos = 1
    seq2_pos = 1

    for i in range(0, len(seq1), line_length):
        seq1_line = seq1[i:i+line_length]
        seq2_line = seq2[i:i+line_length]
        
        # マッチしていない部分を赤色で強調
        seq1_colored = ''.join([f"<span style='color:red'>{base}</span>" if base1 != base2 else base for base1, base2, base in zip(seq1_line, seq2_line, seq1_line)])
        seq2_colored = ''.join([f"<span style='color:red'>{base}</span>" if base1 != base2 else base for base1, base2, base in zip(seq1_line, seq2_line, seq2_line)])

        seq1_end_pos = seq1_pos + len(seq1_line) - 1
        seq2_end_pos = seq2_pos + len(seq2_line) - 1

        formatted_result += f"<pre>read_ {seq1_pos:>4} > {seq1_colored:<{line_length}} < {seq1_end_pos:>4}</pre>"
        formatted_result += f"<pre>target {seq2_pos:>4} > {seq2_colored:<{line_length}} < {seq2_end_pos:>4}</pre>"
        formatted_result += "\n"

        seq1_pos += len(seq1_line)
        seq2_pos += len(seq2_line)

    return formatted_result



