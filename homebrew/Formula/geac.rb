class Geac < Formula
  desc "Genomic Evidence Atlas of Cohorts — collect alt-base metrics and explore coverage"
  homepage "https://github.com/fleharty/GEAC"
  version "FORMULA_VERSION"

  on_macos do
    on_arm do
      url "https://github.com/fleharty/GEAC/releases/download/v#{version}/geac-macos-arm64.tar.gz"
      sha256 "SHA256_MACOS_ARM64"
    end
    on_intel do
      url "https://github.com/fleharty/GEAC/releases/download/v#{version}/geac-macos-x86_64.tar.gz"
      sha256 "SHA256_MACOS_X86_64"
    end
  end

  on_linux do
    on_intel do
      url "https://github.com/fleharty/GEAC/releases/download/v#{version}/geac-linux-x86_64.tar.gz"
      sha256 "SHA256_LINUX_X86_64"
    end
    on_arm do
      url "https://github.com/fleharty/GEAC/releases/download/v#{version}/geac-linux-arm64.tar.gz"
      sha256 "SHA256_LINUX_ARM64"
    end
  end

  depends_on "htslib"
  depends_on "python@3.12"

  resource "geac-apps" do
    url "https://github.com/fleharty/GEAC/archive/refs/tags/v#{version}.tar.gz"
    sha256 "SHA256_SOURCE"
  end

  def install
    bin.install "geac"

    resource("geac-apps").stage do
      libexec.install Dir["app/*"]
    end

    python = Formula["python@3.12"].opt_bin/"python3.12"
    venv = libexec/"venv"
    system python, "-m", "venv", venv
    system venv/"bin/pip", "install", "--upgrade", "pip", "--quiet"
    system venv/"bin/pip", "install", "--quiet",
      "streamlit>=1.35",
      "duckdb>=1.0",
      "altair>=5.3",
      "pandas>=2.0",
      "numpy>=1.26",
      "scipy>=1.13",
      "google-cloud-storage",
      "google-auth"

    (bin/"geac-explorer").write <<~SH
      #!/bin/bash
      exec "#{libexec}/venv/bin/streamlit" run "#{libexec}/geac_explorer.py" "$@"
    SH
    chmod 0755, bin/"geac-explorer"

    (bin/"geac-coverage-explorer").write <<~SH
      #!/bin/bash
      exec "#{libexec}/venv/bin/streamlit" run "#{libexec}/geac_coverage_explorer.py" "$@"
    SH
    chmod 0755, bin/"geac-coverage-explorer"
  end

  test do
    assert_match version.to_s, shell_output("#{bin}/geac --version")
  end
end
