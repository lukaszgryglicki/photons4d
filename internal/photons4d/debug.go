//go:build debug
// +build debug

package photons4d

import (
	"fmt"
	"sync"
)

func DebugLog(format string, args ...interface{}) {
	fmt.Printf("[DEBUG] "+format+"\n", args...)
}

var once sync.Once

func DebugLogOnce(format string, args ...interface{}) {
	once.Do(func() {
		fmt.Printf("[DEBUG] "+format+"\n", args...)
	})
}
