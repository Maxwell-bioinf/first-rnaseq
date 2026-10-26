# Читаем test.tsv: в каждой строке имя гена, таб, число прочтений
counts = {}
with open("test.tsv") as f:
    for line in f:
        gene, n = line.strip().split("\t")
        counts[gene] = int(n)

print("Генов:", len(counts))
best = max(counts, key=counts.get)
print("Максимум:", best, counts[best])

