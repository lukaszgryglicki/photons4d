//go:build !debug
// +build !debug

package photons4d

func DebugLog(format string, args ...interface{}) {
	// no-op in release builds
}

func DebugLogOnce(format string, args ...interface{}) {
	// no-op in release builds
}
