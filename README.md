See [run.sh](run.sh) to execute the code associated with the paper.

## Example of scoring

The following demonstrates how to score the deleteriousness of a mutation using our exchangeability and surface accessibility model.
The working directory here is the root of this repository.
First we make imports, including a function from [structural_features.py](data/structural_features.py) to extract surface accessibilities from a given PDB file:

```python
import math
import pandas as pd
from data.structural_features import extract_accessibilities
```

We then load the parameters of the model:
```python
parameters = pd.read_csv("results/full.csv")
beta_accessibility = parameters.at[0, "sa"]
beta_exchange = {
    (row["wildtype"], row["mutation"]): row["mu"] for _, row in parameters.iterrows()
}
```

And compute the score of the mutation M115H in beta-lactamase:
```python
def score(wildtype_aa: str, mutant_aa: str, accessibility: float) -> float:
    return 1 / (
        1
        + math.exp(
            -beta_exchange[(wildtype_aa, mutant_aa)]
            - accessibility * beta_accessibility
        )
    )

accessibilities = extract_accessibilities("data/structures/BLAT_ECOLX.pdb")
print(score("M", "H", accessibilities[115]))
```

Higher scores signify greater deleteriousness.
