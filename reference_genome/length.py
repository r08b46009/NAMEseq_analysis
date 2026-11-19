from Bio import SeqIO

def calculate_genome_length(fasta_file):
    total_length = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_length += len(record.seq)
    return total_length

# 使用範例
fasta_file = "GCF.fasta"  # 替換成你的 FASTA 檔案名稱
genome_length = calculate_genome_length(fasta_file)
print(f"Total genome length: {genome_length} bp")
