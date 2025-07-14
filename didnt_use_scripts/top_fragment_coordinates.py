from collections import Counter

pairs = Counter(zip(starts, ends))
top_fragments = pairs.most_common(10)

print("Top 10 (start, end) fragment coordinates:")
for (s, e), count in top_fragments:
    print(f"{s}-{e} → {count} reads")
