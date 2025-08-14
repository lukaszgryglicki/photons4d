package photons4d

var (
	Debug           = false // set to true for verbose debug output
	UseLocks        = true  // set to false to disable locks for parallel writes
	PNG             = false // set to true to save PNG sequence of the scene (16-bit per channel, lossless PNGs)
	RAW             = false // set to true to save RAW scene
	AlwaysBVH       = false // set to true to always use BVH for nearest hit calculations
	NeverBVH        = false // set to true to never use BVH for nearest hit calculations
	EscapeSPPAdjust = false // set to true to adjust the scene SPP for environmental 4D sphere scene shape
	NearestHitFunc  = nearestHitBVH
	// Compile time checks to ensure that the material interface is implemented by all required types
	_ material = (*Cell8)(nil)
	_ material = (*HyperSphere)(nil)
	_ material = (*Cell5)(nil)
	_ material = (*Cell16)(nil)
	_ material = (*Cell24)(nil)
	_ material = (*cellPoly)(nil) // for 120/600
)
