import type { Metadata } from "next";
import { Inter, JetBrains_Mono } from "next/font/google";
import { Github } from "lucide-react";
import "./globals.css";

const inter = Inter({
  subsets: ["latin"],
  variable: "--font-inter",
  display: "swap",
});

const jetbrainsMono = JetBrains_Mono({
  subsets: ["latin"],
  variable: "--font-mono-jet",
  display: "swap",
});

export const metadata: Metadata = {
  title: "CADD Pipeline-Novoo",
  description:
    "Compound screening pipeline: analogue generation, ADMET filtering, molecular docking, and retrosynthesis.",
};

const REPO_URL = "https://github.com/Novoo5/cadd-pipeline";

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en" className={`dark ${inter.variable} ${jetbrainsMono.variable}`}>
      <body className="font-sans bg-gray-950 text-gray-200 antialiased min-h-screen flex flex-col">
        {/* ── Nav ──────────────────────────────────────────────── */}
        <header className="sticky top-0 z-50 border-b border-white/[0.07] bg-gray-950/80 backdrop-blur-xl">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
            <div className="flex items-center justify-between h-14">
              <div className="flex items-center gap-2.5">
                {/* Hexagonal mark — a nod to a benzene ring */}
                <svg
                  viewBox="0 0 24 24"
                  className="w-7 h-7 flex-none"
                  aria-hidden="true"
                >
                  <path
                    d="M12 2.6 20.1 7.3v9.4L12 21.4 3.9 16.7V7.3z"
                    fill="rgb(22 212 142 / 0.12)"
                    stroke="currentColor"
                    strokeWidth="1.4"
                    className="text-emerald-500"
                    strokeLinejoin="round"
                  />
                  <circle cx="12" cy="12" r="2.6" className="fill-emerald-400" />
                </svg>
                <div className="flex items-baseline gap-2">
                  <span className="font-semibold text-gray-100 tracking-tight text-sm">
                    CADD Pipeline
                  </span>
                  <span className="hidden sm:inline text-[10px] font-mono text-gray-600 border border-white/[0.07] rounded px-1.5 py-0.5">
                    v1.0
                  </span>
                </div>
              </div>

              <a
                href={REPO_URL}
                target="_blank"
                rel="noopener noreferrer"
                className="inline-flex items-center gap-1.5 text-xs font-medium text-gray-500 hover:text-gray-200 transition-colors"
              >
                <Github className="w-3.5 h-3.5" />
                <span className="hidden sm:inline">Source</span>
              </a>
            </div>
          </div>
        </header>

        {/* ── Content ──────────────────────────────────────────── */}
        <main className="flex-1 max-w-7xl w-full mx-auto px-4 sm:px-6 lg:px-8 py-10 sm:py-14">
          {children}
        </main>

        {/* ── Footer ───────────────────────────────────────────── */}
        <footer className="border-t border-white/[0.07] mt-16">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6">
            <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-3 text-xs text-gray-600">
              <p>
                Built with RDKit, ADMET-AI, DeepChem and AutoDock Vina — all
                open source.
              </p>
              <p className="font-mono text-[11px]">
                FastAPI · Next.js · Hugging Face
              </p>
            </div>
          </div>
        </footer>
      </body>
    </html>
  );
}
