#python

#Подключаем библиотеку numpy и pandas
import numpy as np
import pandas as pd

#Чтение таблицы: колонки разделены табом (\), строки с названиями нет
df = pd.read_csv("test.tsv", sep="\t", header=None, names=["gene", "reads"])

df["cpm"] = df["reads"] / df["reads"].sum() * 1_000_000
df["log_cpm"] = np.log2(df["cpm"] + 1)
print(df)
