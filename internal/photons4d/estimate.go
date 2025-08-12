package photons4d

import (
	"math/rand"
	"runtime"
	"sync"
	"time"
)

func estimateHitProb(light *Light, scene *Scene, trials int) Real {
	if trials <= 0 {
		return 0
	}
	workers := runtime.NumCPU()
	if workers < 1 {
		workers = 1
	}
	if workers > trials {
		workers = trials
	}

	per, rem := trials/workers, trials%workers
	var wg sync.WaitGroup
	hitsCh := make(chan int, workers)

	for w := 0; w < workers; w++ {
		n := per
		if w < rem {
			n++
		}
		if n == 0 {
			continue
		}
		wg.Add(1)
		go func(wid, n int) {
			defer wg.Done()
			// independent RNG per worker
			seed := time.Now().UnixNano() ^ int64(uint64(wid)*0x9e3779b97f4a7c15)
			rng := rand.New(rand.NewSource(seed))

			localHits := 0
			for i := 0; i < n; i++ {
				if castSingleRay(light, scene, rng, nil, false) {
					localHits++
				}
			}
			hitsCh <- localHits
		}(w, n)
	}

	wg.Wait()
	close(hitsCh)

	totalHits := 0
	for h := range hitsCh {
		totalHits += h
	}
	return Real(totalHits) / Real(trials)
}
