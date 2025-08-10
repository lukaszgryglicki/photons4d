//go:build !debug
// +build !debug

package main

func debugLog(format string, args ...interface{}) {
	// no-op in release builds
}

func debugLogOnce(format string, args ...interface{}) {
	// no-op in release builds
}
