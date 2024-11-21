from sub_folder.function import *
import streamlit as st
import os
import re
from io import StringIO
from Bio import SeqIO
from Bio import pairwise2

# Main Streamlitアプリケーション
def main():
    st.title("シーケンスファイルのトリミングとペアワイズアラインメント")

    # GUIからパラメータを取得
    ran5_prime = st.number_input("5'末端のトリミング範囲 (bp):5'末端から指定した領域内で最後にNが出現してから10bp下流までを除きます", value=100)
    read_len = st.number_input("リード長 (bp):リード長より下流でNが出現した部分の上流10bp以降を除きます", value=800)

    # ファイルのアップロード
    uploaded_files = st.file_uploader("シーケンスファイル(read_)を選択してください", accept_multiple_files=True, type=["fasta"])
    comparison_file = st.file_uploader("比較配列(target)", type=["fasta", "gcc"], key="comparison")

    if uploaded_files and comparison_file and st.button('実行'):
        try:
            # アップロードされたファイルを読み込み、トリミングを実行
            reads = []
            for uploaded_file in uploaded_files:
                stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                seq_record = SeqIO.read(stringio, "fasta")
                trimmed_sequence = trim_sequence(str(seq_record.seq), ran5_prime=ran5_prime, read_len=read_len)
                reads.append((uploaded_file.name, trimmed_sequence))

            # 比較配列を読み込む
            comparison_stringio = StringIO(comparison_file.getvalue().decode("utf-8"))
            if comparison_file.name.endswith(".fasta"):
                comparison_seq = str(SeqIO.read(comparison_stringio, "fasta").seq)
            elif comparison_file.name.endswith(".gcc"):
                lines = comparison_stringio.readlines()
                comparison_seq = "".join([line.strip() for line in lines if not line.startswith('>')])
                comparison_seq = re.sub(r'[^ATGC]', '', comparison_seq)

            # 各リードを比較配列にアラインメント
            alignments = [(name, pairwise2.align.globalms(read, comparison_seq, 2, -1, -0.5, -0.1, one_alignment_only=True)[0]) for name, read in reads]

            # 結果をタブで表示
            st.write("Alignment Results:")
            tab_names = [name for name, _ in alignments]
            tabs = st.tabs(tab_names)
            
            for tab, (name, alignment) in zip(tabs, alignments):
                with tab:
                    formatted_result = format_alignment(alignment.seqA, alignment.seqB, read_name=name, target_name="target")
                    st.markdown(formatted_result, unsafe_allow_html=True)

            st.success("トリミングとペアワイズアラインメントが完了しました。")
        
        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    main()
