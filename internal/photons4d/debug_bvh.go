package photons4d

import (
	"fmt"
	"strings"
)

// DumpAABBBVH prints the BVH tree with indentation (one space per level).
// It prints subtree counts (nodes, leaves, objects) and the AABB min/max for each node.
func DumpAABBBVH(scene *Scene, build bool) bool {
	var root *AABBNode
	if build {
		root = getOrBuildBVH(scene)
	} else {
		data, ok := bvhCache.Load(scene)
		if data == nil || !ok {
			fmt.Println("[BVH] <empty>")
			return false
		}
		root = data.(*AABBNode)
	}
	memo := make(map[*AABBNode]bvhCounts, 1024)
	totals := bvhCount(root, memo)
	fmt.Printf("[BVH] root: nodes=%d leaves=%d objs=%d\n", totals.nodes, totals.leaves, totals.objs)
	bvhPrint(root, 0, memo)
	return true
}

type bvhCounts struct {
	nodes  int
	leaves int
	objs   int
}

func bvhCount(n *AABBNode, memo map[*AABBNode]bvhCounts) bvhCounts {
	if n == nil {
		return bvhCounts{}
	}
	if c, ok := memo[n]; ok {
		return c
	}
	if n.leafObjs != nil {
		c := bvhCounts{nodes: 1, leaves: 1, objs: len(n.leafObjs)}
		memo[n] = c
		return c
	}
	lc := bvhCount(n.left, memo)
	rc := bvhCount(n.right, memo)
	c := bvhCounts{
		nodes:  1 + lc.nodes + rc.nodes,
		leaves: lc.leaves + rc.leaves,
		objs:   lc.objs + rc.objs,
	}
	memo[n] = c
	return c
}

func bvhPrint(n *AABBNode, depth int, memo map[*AABBNode]bvhCounts) {
	if n == nil {
		return
	}
	ind := strings.Repeat("\t", depth)
	c := memo[n]
	if n.leafObjs != nil {
		fmt.Printf("%sLEAF  objs=%d | min=(%.5g,%.5g,%.5g,%.5g) max=(%.5g,%.5g,%.5g,%.5g)\n",
			ind, len(n.leafObjs),
			n.min.X, n.min.Y, n.min.Z, n.min.W,
			n.max.X, n.max.Y, n.max.Z, n.max.W,
		)
		return
	}
	fmt.Printf("%sNODE  nodes=%d leaves=%d objs=%d | min=(%.5g,%.5g,%.5g,%.5g) max=(%.5g,%.5g,%.5g,%.5g)\n",
		ind, c.nodes, c.leaves, c.objs,
		n.min.X, n.min.Y, n.min.Z, n.min.W,
		n.max.X, n.max.Y, n.max.Z, n.max.W,
	)
	bvhPrint(n.left, depth+1, memo)
	bvhPrint(n.right, depth+1, memo)
}
