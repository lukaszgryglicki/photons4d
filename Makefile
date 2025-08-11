SRC = photons4d.go
BIN = photons4d

all: build

build: fmt
	go build -tags debug -o $(BIN) $(SRC) debug.go
	strip $(BIN)

release: fmt
	go build -o $(BIN) $(SRC) nodebug.go
	strip -s $(BIN)

run: build
	./$(BIN)

fmt:
	go fmt $(SRC) debug.go nodebug.go

clean:
	rm -f $(BIN)
