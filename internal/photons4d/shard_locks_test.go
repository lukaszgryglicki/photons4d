package photons4d

import "testing"

func TestShardLocksBasic(t *testing.T) {
	sl := &shardLocks{}
	idx := 123456
	sl.lock(idx)
	sl.unlock(idx)
}
