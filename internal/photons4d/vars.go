package photons4d

var (
	Debug          = false // set to true for verbose debug output
	UseLocks       = true  // set to false to disable locks for parallel writes
	PNG            = false // set to true to save PNG sequence of the scene (16-bit per channel, lossless PNGs)
	RAW            = false // set to true to save RAW scene
	NearestHitFunc = nearestHitBVH
)
