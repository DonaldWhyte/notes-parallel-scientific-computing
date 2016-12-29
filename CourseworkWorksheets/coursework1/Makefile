default: build

build:
	$(MAKE) -C parallel build
	$(MAKE) -C parallel_zoom build
	$(MAKE) -C serial build

clean:
	$(MAKE) -C parallel clean
	$(MAKE) -C parallel_zoom clean
	$(MAKE) -C serial clean

.PHONY: clean