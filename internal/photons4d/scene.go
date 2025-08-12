package photons4d

// Scene stores a 3D volume (axis-aligned in X,Y,Z) embedded at W=Center.W.
type Scene struct {
	Center               Point4
	Width, Height, Depth Real
	Nx, Ny, Nz           int
	MaxBounces           int
	Buf                  []Real // flat: (((i*Ny)+j)*Nz + k)*3 + c
	Cells8               []*Cell8
	Hyperspheres         []*HyperSphere
	Cells5               []*Cell5
	Cells16              []*Cell16
	Cells24              []*Cell24
	Cells120             []*Cell120
	Cells600             []*Cell600

	// cached bounds & mapping
	MinX, MaxX Real
	MinY, MaxY Real
	MinZ, MaxZ Real
	InvSpanX   Real
	InvSpanY   Real
	InvSpanZ   Real
	StrideX    int // i * StrideX + j * StrideY + k*3 + c
	StrideY    int
}

// NewScene allocates a zero-initialized flat voxel grid and precomputes bounds & strides.
func NewScene(center Point4, width, height, depth Real, nx, ny, nz, maxBounces int) *Scene {
	if nx <= 0 || ny <= 0 || nz <= 0 {
		panic("voxel resolution must be positive")
	}
	total := nx * ny * nz * 3

	halfX := width * 0.5
	halfY := height * 0.5
	halfZ := depth * 0.5
	minX, maxX := center.X-halfX, center.X+halfX
	minY, maxY := center.Y-halfY, center.Y+halfY
	minZ, maxZ := center.Z-halfZ, center.Z+halfZ

	invSpanX := 1.0 / (maxX - minX)
	invSpanY := 1.0 / (maxY - minY)
	invSpanZ := 1.0 / (maxZ - minZ)

	strideY := nz * 3
	strideX := ny * strideY

	s := &Scene{
		Center:     center,
		Width:      width,
		Height:     height,
		Depth:      depth,
		Nx:         nx,
		Ny:         ny,
		Nz:         nz,
		Buf:        make([]Real, total),
		MaxBounces: maxBounces,

		MinX:     minX,
		MaxX:     maxX,
		MinY:     minY,
		MaxY:     maxY,
		MinZ:     minZ,
		MaxZ:     maxZ,
		InvSpanX: invSpanX,
		InvSpanY: invSpanY,
		InvSpanZ: invSpanZ,
		StrideX:  strideX,
		StrideY:  strideY,
	}
	DebugLog("Created scene center=%+v, size=(%.2f, %.2f, %.2f), resolution=(%d, %d, %d), maxBounces=%d", center, width, height, depth, nx, ny, nz, maxBounces)
	return s
}

// VoxelSize returns the physical size of each voxel along X,Y,Z.
func (s *Scene) VoxelSize() (dx, dy, dz Real) {
	dx, dy, dz = s.Width/Real(s.Nx), s.Height/Real(s.Ny), s.Depth/Real(s.Nz)
	DebugLogOnce("Voxel size: (%.5f, %.5f, %.5f)", dx, dy, dz)
	return
}

// VoxelIndexOf maps a 4D point to voxel indices and also returns normalized coords (u,v,w) in [0,1].
func (s *Scene) VoxelIndexOf(p Point4) (ok bool, i, j, k int, ux, uy, uz Real) {
	if p.X < s.MinX || p.X >= s.MaxX || p.Y < s.MinY || p.Y >= s.MaxY || p.Z < s.MinZ || p.Z >= s.MaxZ {
		return false, 0, 0, 0, 0, 0, 0
	}
	ux = (p.X - s.MinX) * s.InvSpanX
	uy = (p.Y - s.MinY) * s.InvSpanY
	uz = (p.Z - s.MinZ) * s.InvSpanZ
	i = int(ux * Real(s.Nx))
	j = int(uy * Real(s.Ny))
	k = int(uz * Real(s.Nz))
	if i == s.Nx {
		i = s.Nx - 1
	}
	if j == s.Ny {
		j = s.Ny - 1
	}
	if k == s.Nz {
		k = s.Nz - 1
	}
	return true, i, j, k, ux, uy, uz
}

// Flat buffer index helper (c ∈ {ChR,ChG,ChB}).
func (s *Scene) idx(i, j, k, c int) int {
	return i*s.StrideX + j*s.StrideY + k*3 + c
}

func (s *Scene) AddCell8(h *Cell8) {
	s.Cells8 = append(s.Cells8, h)
}

func (s *Scene) AddHyperSphere(h *HyperSphere) {
	s.Hyperspheres = append(s.Hyperspheres, h)
}

func (s *Scene) AddCell5(h *Cell5) {
	s.Cells5 = append(s.Cells5, h)
}

func (s *Scene) AddCell16(h *Cell16) {
	s.Cells16 = append(s.Cells16, h)
}
func (s *Scene) AddCell24(h *Cell24) {
	s.Cells24 = append(s.Cells24, h)
}

func (s *Scene) AddCell120(c *Cell120) {
	s.Cells120 = append(s.Cells120, c)
}
func (s *Scene) AddCell600(c *Cell600) {
	s.Cells600 = append(s.Cells600, c)
}
