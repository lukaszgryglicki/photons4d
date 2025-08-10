//go:build debug
// +build debug

package main

import (
	"fmt"
	"sync"
)

func debugLog(format string, args ...interface{}) {
	fmt.Printf("[DEBUG] "+format+"\n", args...)
}

var once sync.Once

func debugLogOnce(format string, args ...interface{}) {
	once.Do(func() {
		fmt.Printf("[DEBUG] "+format+"\n", args...)
	})
}
