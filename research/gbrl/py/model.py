import torch
import torch.nn as nn


# Pointer-network-style policy: each critical pair is embedded independently
# (concatenated with broadcast global state) and the actor produces one logit
# per pair, so the action distribution is over a variable number of pairs.
# Critic pools the pair embeddings into a single state value.
class PointerPolicy(nn.Module):
    def __init__(self, pair_dim: int = 10, global_dim: int = 3, hidden: int = 64):
        super().__init__()
        in_dim = pair_dim + global_dim
        self.encoder = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.Tanh(),
            nn.Linear(hidden, hidden),
            nn.Tanh(),
        )
        self.logit_head = nn.Linear(hidden, 1)
        self.value_head = nn.Sequential(
            nn.Linear(hidden, hidden),
            nn.Tanh(),
            nn.Linear(hidden, 1),
        )
        # Small init on logit head so initial policy starts near uniform but
        # leaves room for gradients to push pairs apart cleanly.
        nn.init.orthogonal_(self.logit_head.weight, gain=0.01)
        nn.init.zeros_(self.logit_head.bias)
        for m in self.encoder:
            if isinstance(m, nn.Linear):
                nn.init.orthogonal_(m.weight, gain=1.0)
                nn.init.zeros_(m.bias)
        for m in self.value_head:
            if isinstance(m, nn.Linear):
                nn.init.orthogonal_(m.weight, gain=1.0)
                nn.init.zeros_(m.bias)

    def forward(self, pair_feats: torch.Tensor, global_feats: torch.Tensor):
        # Single-sample forward. pair_feats: [n_pairs, P], global_feats: [G].
        n = pair_feats.shape[0]
        g = global_feats.unsqueeze(0).expand(n, -1)
        x = torch.cat([pair_feats, g], dim=-1)
        emb = self.encoder(x)                      # [n_pairs, H]
        logits = self.logit_head(emb).squeeze(-1)  # [n_pairs]
        pooled = emb.mean(dim=0)                   # [H]
        value = self.value_head(pooled).squeeze(-1)  # scalar
        return logits, value

    def forward_batched(
        self,
        pair_feats: torch.Tensor,    # [B, P_max, P]
        global_feats: torch.Tensor,  # [B, G]
        pair_mask: torch.Tensor,     # [B, P_max] bool, True for valid pairs
    ):
        """Batched forward. Pad pair_feats to P_max along dim 1; pair_mask
        marks valid positions. Logits at padded positions are set to -inf so
        Categorical(logits=...) ignores them. Value pools over valid pairs."""
        B, P_max, _ = pair_feats.shape
        g = global_feats.unsqueeze(1).expand(B, P_max, -1)
        x = torch.cat([pair_feats, g], dim=-1)
        emb = self.encoder(x)                                # [B, P_max, H]
        logits = self.logit_head(emb).squeeze(-1)            # [B, P_max]
        logits = logits.masked_fill(~pair_mask, float('-inf'))
        # Mean over valid pairs (avoid div-by-zero with .clamp).
        mask_f = pair_mask.float().unsqueeze(-1)             # [B, P_max, 1]
        pooled = (emb * mask_f).sum(dim=1) / mask_f.sum(dim=1).clamp(min=1.0)
        value = self.value_head(pooled).squeeze(-1)          # [B]
        return logits, value
