package main

import (
	"fmt"
	"math/rand"
	"os"
	"os/signal"
	"runtime/pprof"
	"sync"
	"syscall"
	"time"

	"github.com/lukaszgryglicki/photons4d/internal/photons4d"
)

func main() {
	rand.Seed(time.Now().UnixNano())
	cleanup := func() {}
	var cleanupOnce sync.Once
	runCleanup := func() {
		cleanupOnce.Do(cleanup)
	}
	defer runCleanup()

	sigCh := make(chan os.Signal, 1)
	signal.Notify(sigCh, syscall.SIGUSR1)
	defer signal.Stop(sigCh)
	go func() {
		<-sigCh
		fmt.Fprintln(os.Stderr, "[INFO] SIGUSR1 received, exiting successfully")
		runCleanup()
		os.Exit(0)
	}()

	photons4d.Debug = os.Getenv("DEBUG") != ""
	photons4d.UseLocks = os.Getenv("SKIP_LOCKS") == ""
	photons4d.PNG = os.Getenv("PNG") != ""
	photons4d.RAW = os.Getenv("RAW") != ""
	photons4d.AlwaysBVH = os.Getenv("ALWAYS_BVH") != ""
	photons4d.NeverBVH = os.Getenv("NEVER_BVH") != ""
	photons4d.EscapeSPPAdjust = os.Getenv("SPP_ADJUST") != ""
	photons4d.ForceEscape = os.Getenv("FORCE_ESCAPE") != ""
	profile := os.Getenv("PROFILE") != ""
	if profile {
		f, err := os.Create("cpu.out")
		if err != nil {
			panic(err)
		}
		if err := pprof.StartCPUProfile(f); err != nil {
			panic(err)
		}
		cleanup = func() {
			pprof.StopCPUProfile()
			_ = f.Close()
		}
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
