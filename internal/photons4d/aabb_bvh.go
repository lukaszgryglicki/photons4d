package photons4d

import (
	"math"
	"sort"
	"sync"
)

type bvhLeaf struct {
	min, max  Point4
	intersect func(Point4, Vector4) (objectHit, bool)
}

type AABBNode struct {
	min, max Point4
	left     *AABBNode
	right    *AABBNode
	leafObjs []bvhLeaf // non-nil ⇒ leaf
}

// BVH cache keyed by *Scene (avoids changing Scene struct).
var bvhCache sync.Map // map[*Scene]*AABBNode

func getOrBuildBVH(s *Scene) *AABBNode {
	if v, ok := bvhCache.Load(s); ok {
		return v.(*AABBNode)
	}
	objs := collectSceneObjects(s)
	var root *AABBNode
	if len(objs) > 0 {
		root = buildBVH(objs)
	} else {
		// Empty scene ⇒ nil root (no geometry)
		root = nil
	}
	bvhCache.Store(s, root)
	return root
}

func collectSceneObjects(s *Scene) []bvhLeaf {
	out := make([]bvhLeaf, 0,
		len(s.Cells8)+len(s.Hyperspheres)+len(s.Cells5)+len(s.Cells16)+len(s.Cells24)+len(s.Cells120)+len(s.Cells600))

	for _, o := range s.Cells8 {
		if o == nil {
			continue
		}
		obj := o
		out = append(out, bvhLeaf{
			min: obj.AABBMin, max: obj.AABBMax,
			intersect: func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell8(O, D, obj) },
		})
	}
	for _, o := range s.Hyperspheres {
		if o == nil {
			continue
		}
		obj := o
		out = append(out, bvhLeaf{
			min: obj.AABBMin, max: obj.AABBMax,
			intersect: func(O Point4, D Vector4) (objectHit, bool) { return intersectRayHyperSphere(O, D, obj) },
		})
	}
	for _, o := range s.Cells5 {
		if o == nil {
			continue
		}
		obj := o
		out = append(out, bvhLeaf{
			min: obj.AABBMin, max: obj.AABBMax,
			intersect: func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell5(O, D, obj) },
		})
	}
	for _, o := range s.Cells16 {
		if o == nil {
			continue
		}
		obj := o
		out = append(out, bvhLeaf{
			min: obj.AABBMin, max: obj.AABBMax,
			intersect: func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell16(O, D, obj) },
		})
	}
	for _, o := range s.Cells24 {
		if o == nil {
			continue
		}
		obj := o
		out = append(out, bvhLeaf{
			min: obj.AABBMin, max: obj.AABBMax,
			intersect: func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCell24(O, D, obj) },
		})
	}
	for _, o := range s.Cells120 {
		if o == nil {
			continue
		}
		obj := o
		out = append(out, bvhLeaf{
			min: obj.AABBMin, max: obj.AABBMax,
			intersect: func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCellPoly(O, D, &obj.cellPoly) },
		})
	}
	for _, o := range s.Cells600 {
		if o == nil {
			continue
		}
		obj := o
		out = append(out, bvhLeaf{
			min: obj.AABBMin, max: obj.AABBMax,
			intersect: func(O Point4, D Vector4) (objectHit, bool) { return intersectRayCellPoly(O, D, &obj.cellPoly) },
		})
	}
	return out
}

func buildBVH(objs []bvhLeaf) *AABBNode {
	return buildBVHRec(objs, 0)
}

func buildBVHRec(objs []bvhLeaf, depth int) *AABBNode {
	n := len(objs)
	if n == 0 {
		return nil
	}
	if n <= AABBBVHMaxLeafSize {
		minP, maxP := objs[0].min, objs[0].max
		for i := 1; i < n; i++ {
			minP, maxP = aabbUnion(minP, maxP, objs[i].min, objs[i].max)
		}
		return &AABBNode{min: minP, max: maxP, leafObjs: objs}
	}

	// Union bounds and centroid spreads
	minP, maxP := objs[0].min, objs[0].max
	cmin := [4]Real{centroid(objs[0].min.X, objs[0].max.X), centroid(objs[0].min.Y, objs[0].max.Y), centroid(objs[0].min.Z, objs[0].max.Z), centroid(objs[0].min.W, objs[0].max.W)}
	cmax := cmin
	for i := 1; i < n; i++ {
		minP, maxP = aabbUnion(minP, maxP, objs[i].min, objs[i].max)
		cx := centroid(objs[i].min.X, objs[i].max.X)
		cy := centroid(objs[i].min.Y, objs[i].max.Y)
		cz := centroid(objs[i].min.Z, objs[i].max.Z)
		cw := centroid(objs[i].min.W, objs[i].max.W)
		if cx < cmin[0] {
			cmin[0] = cx
		}
		if cx > cmax[0] {
			cmax[0] = cx
		}
		if cy < cmin[1] {
			cmin[1] = cy
		}
		if cy > cmax[1] {
			cmax[1] = cy
		}
		if cz < cmin[2] {
			cmin[2] = cz
		}
		if cz > cmax[2] {
			cmax[2] = cz
		}
		if cw < cmin[3] {
			cmin[3] = cw
		}
		if cw > cmax[3] {
			cmax[3] = cw
		}
	}
	spread := [4]Real{
		cmax[0] - cmin[0],
		cmax[1] - cmin[1],
		cmax[2] - cmin[2],
		cmax[3] - cmin[3],
	}
	axis := 0
	if spread[1] > spread[axis] {
		axis = 1
	}
	if spread[2] > spread[axis] {
		axis = 2
	}
	if spread[3] > spread[axis] {
		axis = 3
	}

	// If all centroids coincide (degenerate), fall back to longest box extent axis.
	if spread[axis] <= 1e-18 {
		ext := [4]Real{
			maxP.X - minP.X,
			maxP.Y - minP.Y,
			maxP.Z - minP.Z,
			maxP.W - minP.W,
		}
		axis = 0
		if ext[1] > ext[axis] {
			axis = 1
		}
		if ext[2] > ext[axis] {
			axis = 2
		}
		if ext[3] > ext[axis] {
			axis = 3
		}
	}

	// Sort by chosen centroid axis, split at median
	sort.Slice(objs, func(i, j int) bool {
		ci := getCentroidAxis(objs[i], axis)
		cj := getCentroidAxis(objs[j], axis)
		if ci == cj {
			return i < j
		}
		return ci < cj
	})
	mid := n / 2
	left := buildBVHRec(objs[:mid], depth+1)
	right := buildBVHRec(objs[mid:], depth+1)

	node := &AABBNode{min: minP, max: maxP, left: left, right: right}
	return node
}

func aabbUnion(aMin, aMax, bMin, bMax Point4) (Point4, Point4) {
	return Point4{
			rmin(aMin.X, bMin.X),
			rmin(aMin.Y, bMin.Y),
			rmin(aMin.Z, bMin.Z),
			rmin(aMin.W, bMin.W),
		}, Point4{
			rmax(aMax.X, bMax.X),
			rmax(aMax.Y, bMax.Y),
			rmax(aMax.Z, bMax.Z),
			rmax(aMax.W, bMax.W),
		}
}

func centroid(a, b Real) Real { return (a + b) * 0.5 }

func getCentroidAxis(o bvhLeaf, axis int) Real {
	switch axis {
	case 0:
		return centroid(o.min.X, o.max.X)
	case 1:
		return centroid(o.min.Y, o.max.Y)
	case 2:
		return centroid(o.min.Z, o.max.Z)
	default:
		return centroid(o.min.W, o.max.W)
	}
}

func computeRayRecips(d Vector4) rayRecips {
	const eps = 1e-18
	rr := rayRecips{}
	if x := d.X; x > eps || x < -eps {
		rr.invX = 1 / x
	} else {
		rr.parX = true
	}
	if y := d.Y; y > eps || y < -eps {
		rr.invY = 1 / y
	} else {
		rr.parY = true
	}
	if z := d.Z; z > eps || z < -eps {
		rr.invZ = 1 / z
	} else {
		rr.parZ = true
	}
	if w := d.W; w > eps || w < -eps {
		rr.invW = 1 / w
	} else {
		rr.parW = true
	}
	return rr
}

// Nearest-hit traversal (iterative, stack-based). Prunes by current best t.
func traverseNearest(root *AABBNode, O Point4, D Vector4, tMax Real) (objectHit, bool) {
	if root == nil {
		return objectHit{}, false
	}
	bestT := tMax
	var best objectHit
	rr := computeRayRecips(D)

	type entry struct {
		n    *AABBNode
		tmin Real
	}
	stack := []entry{{n: root, tmin: 0}}
	for len(stack) > 0 {
		// pop
		e := stack[len(stack)-1]
		stack = stack[:len(stack)-1]

		ok, tmin := rayAABB(O, e.n.min, e.n.max, rr)
		if !ok || tmin > bestT || tmin < 0 {
			continue
		}

		if e.n.leafObjs != nil {
			for i := range e.n.leafObjs {
				if h, ok := e.n.leafObjs[i].intersect(O, D); ok && h.t > 0 && h.t < bestT {
					bestT = h.t
					best = h
				}
			}
			continue
		}

		// order children near→far (push far first so near is processed next)
		var lOK bool
		var lT Real
		if e.n.left != nil {
			lOK, lT = rayAABB(O, e.n.left.min, e.n.left.max, rr)
			lOK = lOK && lT <= bestT
		}
		var rOK bool
		var rT Real
		if e.n.right != nil {
			rOK, rT = rayAABB(O, e.n.right.min, e.n.right.max, rr)
			rOK = rOK && rT <= bestT
		}
		if lOK && rOK {
			if lT < rT {
				stack = append(stack, entry{e.n.right, rT}, entry{e.n.left, lT})
			} else {
				stack = append(stack, entry{e.n.left, lT}, entry{e.n.right, rT})
			}
		} else if lOK {
			stack = append(stack, entry{e.n.left, lT})
		} else if rOK {
			stack = append(stack, entry{e.n.right, rT})
		}
	}

	if bestT < tMax && bestT < Real(math.Inf(1)) {
		return best, true
	}
	return objectHit{}, false
}
