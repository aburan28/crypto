// Render the complete SVG viewport with native WebKit, preserving its aspect ratio.
import AppKit
import WebKit
let args = CommandLine.arguments
guard args.count == 5, let width = Double(args[3]), let height = Double(args[4]) else { fatalError("SVG PNG WIDTH HEIGHT") }
final class Renderer: NSObject, WKNavigationDelegate {
    let view: WKWebView
    let window: NSWindow
    let output: URL
    init(width: Double, height: Double, output: URL) {
        self.output = output
        view = WKWebView(frame: NSRect(x: 0, y: 0, width: width, height: height))
        window = NSWindow(contentRect: view.frame, styleMask: .borderless, backing: .buffered, defer: false)
        super.init()
        window.contentView = view
        view.navigationDelegate = self
    }
    func webView(_ webView: WKWebView, didFinish navigation: WKNavigation!) {
        DispatchQueue.main.asyncAfter(deadline: .now() + 0.3) {
            let config = WKSnapshotConfiguration()
            config.rect = self.view.bounds
            self.view.takeSnapshot(with: config) { image, error in
                guard let image = image, let tiff = image.tiffRepresentation,
                      let bitmap = NSBitmapImageRep(data: tiff),
                      let png = bitmap.representation(using: .png, properties: [:]) else { fatalError("snapshot failed: \(String(describing: error))") }
                do { try png.write(to: self.output) } catch { fatalError("write failed: \(error)") }
                NSApplication.shared.stop(nil)
                NSApplication.shared.terminate(nil)
            }
        }
    }
}
let app = NSApplication.shared
app.setActivationPolicy(.prohibited)
let renderer = Renderer(width: width, height: height, output: URL(fileURLWithPath: args[2]))
let source = try String(contentsOfFile: args[1], encoding: .utf8)
renderer.view.loadHTMLString("<!doctype html><html><head><meta charset='utf-8'><style>html,body{margin:0;padding:0;background:white;overflow:hidden}</style></head><body>\(source)</body></html>", baseURL: nil)
DispatchQueue.main.asyncAfter(deadline: .now() + 20) { fatalError("render timeout") }
app.run()
