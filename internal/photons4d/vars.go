package photons4d

var (
	Debug          = false // set to true for verbose debug output
	UseLocks       = true  // set to false to disable locks for parallel writes
	PNG            = false // set to true to save PNG sequence of the scene (16-bit per channel, lossless PNGs)
	RAW            = false // set to true to save RAW scene
	NearestHitFunc = nearestHitBVH
	// Compile time checks to ensure that the material interface is implemented by all required types
	_ material = (*Cell8)(nil)
	_ material = (*HyperSphere)(nil)
	_ material = (*Cell5)(nil)
	_ material = (*Cell16)(nil)
	_ material = (*Cell24)(nil)
	_ material = (*cellPoly)(nil) // for 120/600
)
