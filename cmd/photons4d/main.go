package main

import (
	"fmt"
	"math/rand"
	"os"
	"runtime/pprof"
	"time"

	"github.com/lukaszgryglicki/photons4d/internal/photons4d"
)

func main() {
	rand.Seed(time.Now().UnixNano())

	photons4d.Debug = os.Getenv("DEBUG") != ""
	photons4d.UseLocks = os.Getenv("SKIP_LOCKS") == ""
	photons4d.PNG = os.Getenv("PNG") != ""
	photons4d.RAW = os.Getenv("RAW") != ""
	photons4d.AlwaysBVH = os.Getenv("ALWAYS_BVH") != ""
	photons4d.NeverBVH = os.Getenv("NEVER_BVH") != ""
	profile := os.Getenv("PROFILE") != ""
	if profile {
		f, err := os.Create("cpu.out")
		if err != nil {
			panic(err)
		}
		if err := pprof.StartCPUProfile(f); err != nil {
			panic(err)
		}
		defer func() {
			pprof.StopCPUProfile()
			_ = f.Close()
		}()
	}

	cfg := "scenes/config.json"
	if len(os.Args) > 1 {
		cfg = os.Args[1]
	}
	if err := photons4d.Run(cfg); err != nil {
		fmt.Printf("Error: %v\n", err)
		os.Exit(1)
	}
}
