# python

import pandas as pd # pandas 3.x в окружении bioinf

#Чтение таблицы: колонки разделены табом (\), строки с названиями нет
df = pd.read_csv("test.tsv", sep="\t", header=None, names=["gene", "reads"])

print("Генов:", len(df))
best = df.sort_values("reads", ascending=False).iloc[0]
print("Максимум:", best["gene"], best["reads"])

