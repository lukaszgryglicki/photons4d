GO ?= go
BIN = photons4d
CMD = ./cmd/photons4d
INSTALL_DIR ?= /data/scripts
GO_LDFLAGS ?= -s -w
GOAMD64 ?=

# The x86-64 level of an optimized build: the one given as GOAMD64, else the
# highest one that runs on this machine (native amd64 builds only; v4 needs
# AVX-512 and v3 AVX2, and a binary of a higher level does not start at all).
AMD64_LEVEL = level="$(GOAMD64)"; \
	if [ -z "$$level" ] && [ "$$($(GO) env GOARCH)" = amd64 ] && \
		[ "$$($(GO) env GOOS)/$$($(GO) env GOARCH)" = "$$($(GO) env GOHOSTOS)/$$($(GO) env GOHOSTARCH)" ]; then \
		probe="$$(mktemp -d)"; \
		printf 'package main\n\nfunc main() {}\n' > "$$probe/main.go"; \
		for l in v4 v3 v2; do \
			if GOAMD64=$$l $(GO) run "$$probe/main.go" >/dev/null 2>&1; then level=$$l; break; fi; \
		done; \
		rm -rf "$$probe"; \
	fi; \
	if [ -n "$$level" ]; then GOAMD64=$$level; export GOAMD64; echo "GOAMD64=$$level"; fi;

.PHONY: all build release static install run fmt race runprofile profile getprofile test clean

all: build

build: fmt
	$(GO) build -tags debug -o $(BIN) $(CMD)

release: fmt
	@$(AMD64_LEVEL) set -x; $(GO) build -ldflags "$(GO_LDFLAGS)" -o $(BIN) $(CMD)

# A pure Go build without cgo is a static executable on Linux and FreeBSD.
# macOS has no static executables: there it must need nothing but the
# libraries of the system.
static: fmt
	@$(AMD64_LEVEL) set -x; CGO_ENABLED=0 $(GO) build -ldflags "$(GO_LDFLAGS)" -o $(BIN).static $(CMD)
	@set -e; \
	case "$$(uname -s)" in \
		Darwin) \
			others="$$(otool -L $(BIN).static | awk 'NR > 1 && $$1 !~ /^\/(usr\/lib|System\/Library)\// { print $$1 }' | tr '\n' ' ')"; \
			[ -z "$$others" ] || { echo "$(BIN).static needs more than the libraries of macOS: $$others" >&2; exit 1; }; \
			echo "$(BIN).static: only the libraries of macOS are dynamic" ;; \
		*) \
			kind="$$(LC_ALL=C file -b $(BIN).static)"; \
			case "$$kind" in \
				*"statically linked"*|*"static-pie linked"*) echo "$(BIN).static: statically linked" ;; \
				*) echo "$(BIN).static is not a static executable: $$kind" >&2; exit 1 ;; \
			esac ;; \
	esac

install: release static
	mkdir -p "$(INSTALL_DIR)"
	install -m 755 $(BIN).static "$(INSTALL_DIR)/$(BIN)"

run: build
	./$(BIN)

runprofile: release
	PROFILE=1 ./$(BIN)

fmt:
	$(GO) fmt ./...

race:
	$(GO) build -race -tags debug -o $(BIN) $(CMD)

profile: runprofile getprofile

getprofile:
	$(GO) tool pprof -text cpu.out > cpu.out.txt && cat cpu.out.txt

test: fmt
	$(GO) test ./internal/photons4d -v

clean:
	rm -f $(BIN) $(BIN).static cpu.out.txt cpu.out out
