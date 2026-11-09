"""Maximum conjugation: choose which bonds are double so that as many atoms as
possible take part in exactly one double bond (a maximum matching, found with
Edmonds' blossom algorithm; no extra dependency).

Pure graph code: ``edit.py`` decides which atoms may take part and applies the
result.
"""

from __future__ import annotations

from collections import deque


def max_matching(n: int, edges: list[tuple[int, int]], initial: list[tuple[int, int]] = ()) -> list[int]:
    """A maximum-cardinality matching of the graph with ``n`` vertices.

    ``initial`` is a valid matching to start from: augmenting never unmatches a
    vertex, so edges of the start matching stay matched wherever possible
    (the result changes as little as it can). Returns ``mate`` with
    ``mate[v]`` the partner of ``v`` or -1.
    """
    adj: list[list[int]] = [[] for _ in range(n)]
    for a, b in edges:
        if a != b and b not in adj[a]:
            adj[a].append(b)
            adj[b].append(a)
    mate = [-1] * n
    for a, b in initial:
        if mate[a] == -1 and mate[b] == -1 and b in adj[a]:
            mate[a], mate[b] = b, a

    def find_path(root: int) -> tuple[int, list[int]]:
        """Search for an augmenting path from ``root``: (free end or -1, parent array)."""
        parent = [-1] * n
        base = list(range(n))
        in_tree = [False] * n
        in_tree[root] = True
        queue = deque([root])

        def lca(a: int, b: int) -> int:
            seen = [False] * n
            while True:
                a = base[a]
                seen[a] = True
                if mate[a] == -1:
                    break
                a = parent[mate[a]]
            while True:
                b = base[b]
                if seen[b]:
                    return b
                b = parent[mate[b]]

        def mark_path(v: int, b: int, child: int, blossom: list[bool]) -> None:
            while base[v] != b:
                blossom[base[v]] = blossom[base[mate[v]]] = True
                parent[v] = child
                child = mate[v]
                v = parent[mate[v]]

        while queue:
            v = queue.popleft()
            for to in adj[v]:
                if base[v] == base[to] or mate[v] == to:
                    continue
                if to == root or (mate[to] != -1 and parent[mate[to]] != -1):
                    # An odd cycle: contract the blossom.
                    current = lca(v, to)
                    blossom = [False] * n
                    mark_path(v, current, to, blossom)
                    mark_path(to, current, v, blossom)
                    for i in range(n):
                        if blossom[base[i]]:
                            base[i] = current
                            if not in_tree[i]:
                                in_tree[i] = True
                                queue.append(i)
                elif parent[to] == -1:
                    parent[to] = v
                    if mate[to] == -1:
                        return to, parent
                    in_tree[mate[to]] = True
                    queue.append(mate[to])
        return -1, parent

    for root in range(n):
        if mate[root] != -1 or not adj[root]:
            continue
        end, parent = find_path(root)
        while end != -1:  # flip the matching along the path
            previous = parent[end]
            following = mate[previous]
            mate[end], mate[previous] = previous, end
            end = following
    return mate
