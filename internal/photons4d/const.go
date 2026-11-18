package photons4d

// Channel indices for readability.
const (
	ChR        = 0
	ChG        = 1
	ChB        = 2
	SceneResX  = 256
	SceneResY  = 256
	SceneResZ  = 256
	ProbeRays  = 100_000
	Spp        = 128 // samples per voxel target
	GIFOut     = "volume.gif"
	GIFDelay   = 5 // 100ths of a second per frame
	Gamma      = 0.75
	MaxBounces = 32
	// NumShards           = 1024
	NumShards           = 65536
	AttenuateD2         = false // if true, attenuate light by 1/d^2 (distance squared) for each bounce (not needed in radiance model)
	CapRejectThreshold  = 0.1   // cone lights covering at least this part of S^3 reject whole-sphere samples, narrower ones sample the polar angle
	AABBBVHMaxLeafSize  = 6
	AABBBVHFromNObjects = 14 // minimum number of objects to use BVH of AABBs, otherwise just iterate all objects on the scene
	EscapeSppFactor     = 4  // factor to increase spp for escape rays to avoid noise
	// hot-loop constants reused across bounces
	epsDist   = 1e-6
	bumpShift = 1e-6
)
