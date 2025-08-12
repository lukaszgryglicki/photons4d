package photons4d

import "sync"

type shardLocks struct{ mu [NumShards]sync.Mutex }

func (sl *shardLocks) lock(idx int)   { sl.mu[idx&(NumShards-1)].Lock() }
func (sl *shardLocks) unlock(idx int) { sl.mu[idx&(NumShards-1)].Unlock() }
