# Diagram sources

These Graphviz sources generate the PNG figures embedded by the parent Isogeny Labs article.

Render all figures locally with:

```bash
mkdir -p ../assets
for f in *.dot; do
  dot -Tpng -Gdpi=170 "$f" -o "../assets/$(basename "$f" .dot).png"
done
```

The source files are intentionally kept beside the published assets so the diagrams remain editable and reproducible.
