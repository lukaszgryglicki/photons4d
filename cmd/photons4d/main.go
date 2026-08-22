package main

import (
	"flag"
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

func envOr(key, def string) string {
	if v := os.Getenv(key); v != "" {
		return v
	}
	return def
}

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

	// Distributed (horizontal) mode is optional; default "local" keeps the
	// classic single-machine vertical behavior byte-for-byte identical.
	mode := flag.String("mode", envOr("MODE", "local"), "run mode: local (default), server (distributed collector) or client (distributed worker)")
	listen := flag.String("listen", envOr("LISTEN", ":31417"), "server mode: TCP address to listen on")
	connect := flag.String("connect", envOr("CONNECT", ""), "client mode: server address host:port")
	chunks := flag.Int("chunks", 64, "server mode: number of work chunks each light's ray budget is split into (or CHUNKS env)")
	round := flag.Float64("round", 0, "server mode: optional adaptive round sizing — target seconds per client round; 0 (default) = classic fixed -chunks distribution (or ROUND env)")
	stats := flag.Float64("stats", 300, "server mode: progress/per-client stats print interval in seconds (or STATS env)")
	compress := flag.Bool("compress", envOr("COMPRESS", "") != "", "client mode: losslessly compress voxel updates (bit-exact, saves network traffic; or COMPRESS env)")
	calgo := flag.String("calgo", envOr("CALGO", "zstd"), "client mode: compression algorithm when -compress is set: zstd (default) or flate (or CALGO env)")
	clevel := flag.Int("clevel", 6, "client mode: compression level 1..9, higher = smaller wire size, more CPU (or CLEVEL env)")
	batch := flag.Int("batch", 0, "client mode: max sparse entries per update message, 0 = default 4194304 ≈ 64 MB raw (or BATCH env)")
	flag.Parse()
	if v := os.Getenv("CHUNKS"); v != "" {
		var c int
		if _, err := fmt.Sscanf(v, "%d", &c); err == nil && c > 0 {
			*chunks = c
		}
	}
	if v := os.Getenv("ROUND"); v != "" {
		var r float64
		if _, err := fmt.Sscanf(v, "%g", &r); err == nil && r >= 0 {
			*round = r
		}
	}
	if v := os.Getenv("BATCH"); v != "" {
		var b int
		if _, err := fmt.Sscanf(v, "%d", &b); err == nil && b > 0 {
			*batch = b
		}
	}
	if v := os.Getenv("CLEVEL"); v != "" {
		var l int
		if _, err := fmt.Sscanf(v, "%d", &l); err == nil && l >= 1 && l <= 9 {
			*clevel = l
		}
	}
	if v := os.Getenv("STATS"); v != "" {
		var s float64
		if _, err := fmt.Sscanf(v, "%f", &s); err == nil && s > 0 {
			*stats = s
		}
	}

	cfg := "scenes/config.json"
	if flag.NArg() > 0 {
		cfg = flag.Arg(0)
	}

	var err error
	switch *mode {
	case "local":
		err = photons4d.Run(cfg)
	case "server":
		err = photons4d.RunServer(cfg, *listen, *chunks, *round, *stats)
	case "client":
		if *connect == "" {
			err = fmt.Errorf("client mode requires -connect host:port (or CONNECT env)")
		} else {
			var codec uint8 = photons4d.CodecZstd
			switch *calgo {
			case "zstd":
			case "flate", "deflate":
				codec = photons4d.CodecFlate
			default:
				fmt.Printf("Error: unknown -calgo %q (want zstd or flate)\n", *calgo)
				os.Exit(1)
			}
			err = photons4d.RunClient(cfg, *connect, photons4d.ClientOpts{Compress: *compress, Codec: codec, Level: *clevel, BatchEntries: *batch})
		}
	default:
		err = fmt.Errorf("unknown mode %q (want local, server or client)", *mode)
	}
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		os.Exit(1)
	}
}
