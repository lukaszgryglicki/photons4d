BIN := photons4d
CMD := ./cmd/photons4d

.PHONY: all build release run fmt runprofile profile clean

all: build

build: fmt
	go build -tags debug -o $(BIN) $(CMD)

release: fmt
	go build -o $(BIN) $(CMD)
	strip $(BIN)

run: build
	./$(BIN)

runprofile: release
	PROFILE=1 ./$(BIN)

fmt:
	go fmt ./...

profile: runprofile
	go tool pprof -text cpu.out > cpu.out.txt && cat cpu.out.txt

clean:
	rm -f $(BIN) cpu.out out
