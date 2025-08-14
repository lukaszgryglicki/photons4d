package photons4d

// Channel indices for readability.
const (
	ChR                 = 0
	ChG                 = 1
	ChB                 = 2
	SceneResX           = 256
	SceneResY           = 256
	SceneResZ           = 256
	ProbeRays           = 100_000
	Spp                 = 128 // samples per voxel target
	GIFOut              = "volume.gif"
	GIFDelay            = 5 // 100ths of a second per frame
	Gamma               = 0.75
	MaxBounces          = 32
	NumShards           = 1024
	AttenuateD2         = false // if true, attenuate light by 1/d^2 (distance squared) for each bounce (not needed in radiance model)
	LUTRejectThreshold  = 0.1   // threshold for using rejection sampling vs. LUT inverse-CDF sampling for cone light S^3 cap sampling
	LutN                = 256
	AABBBVHMaxLeafSize  = 2
	AABBBVHFromNObjects = 8 // minimum number of objects to use BVH of AABBs, otherwise just iterate all objects on the scene
	// hot-loop constants reused across bounces
	epsDist   = 1e-6
	bumpShift = 1e-6
)
