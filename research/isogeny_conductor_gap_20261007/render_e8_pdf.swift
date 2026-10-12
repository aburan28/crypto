import AppKit
import CoreText
import Foundation

guard CommandLine.arguments.count == 2 else { fatalError("expected output PDF path") }
let output = URL(fileURLWithPath: CommandLine.arguments[1])
var pageBox = CGRect(x: 0, y: 0, width: 612, height: 792)
guard let context = CGContext(output as CFURL, mediaBox: &pageBox, nil) else {
    fatalError("could not open PDF context")
}

let ink = NSColor(calibratedRed: 0.10, green: 0.15, blue: 0.22, alpha: 1)
let accent = NSColor(calibratedRed: 0.07, green: 0.39, blue: 0.55, alpha: 1)
let pale = NSColor(calibratedRed: 0.91, green: 0.96, blue: 0.98, alpha: 1)

func text(_ value: String, x: CGFloat, top: CGFloat, width: CGFloat,
          size: CGFloat = 10.5, bold: Bool = false, color: NSColor = ink) -> CGFloat {
    let font = bold ? NSFont.boldSystemFont(ofSize: size) : NSFont.systemFont(ofSize: size)
    let style = NSMutableParagraphStyle()
    style.lineSpacing = 2.2
    let attributed = NSAttributedString(string: value, attributes: [
        .font: font, .foregroundColor: color, .paragraphStyle: style
    ])
    let setter = CTFramesetterCreateWithAttributedString(attributed)
    let suggested = CTFramesetterSuggestFrameSizeWithConstraints(
        setter, CFRangeMake(0, attributed.length), nil,
        CGSize(width: width, height: 10000), nil
    )
    let height = ceil(suggested.height) + 3
    let frame = CTFramesetterCreateFrame(setter, CFRangeMake(0, attributed.length),
        CGPath(rect: CGRect(x: x, y: top - height, width: width, height: height),
               transform: nil), nil)
    context.textMatrix = .identity
    CTFrameDraw(frame, context)
    return top - height
}

func line(_ y: CGFloat) {
    context.setStrokeColor(accent.cgColor)
    context.setLineWidth(0.7)
    context.move(to: CGPoint(x: 45, y: y))
    context.addLine(to: CGPoint(x: 567, y: y))
    context.strokePath()
}

func box(_ rect: CGRect, title: String, body: String) {
    context.setFillColor(pale.cgColor)
    context.fill(rect)
    context.setStrokeColor(accent.cgColor)
    context.stroke(rect)
    _ = text(title, x: rect.minX + 12, top: rect.maxY - 10,
             width: rect.width - 24, size: 11, bold: true, color: accent)
    _ = text(body, x: rect.minX + 12, top: rect.maxY - 36,
             width: rect.width - 24, size: 10)
}

func pageFooter(_ number: Int) {
    line(49)
    _ = text("E8 source and arithmetic preflight  |  2026-10-09", x: 45, top: 40,
             width: 440, size: 8.5, color: accent)
    _ = text("\(number) / 2", x: 520, top: 40, width: 47, size: 8.5, color: accent)
}

context.beginPDFPage(nil)
var y: CGFloat = 745
y = text("E8: characteristic-two Kani preflight", x: 45, top: y,
         width: 522, size: 21, bold: true, color: accent) - 5
y = text("C37 inert conductor step  |  source audit and exact arithmetic",
         x: 45, top: y, width: 522, size: 11.5) - 10
line(y)
y -= 17
y = text("Question. Does Galbraith's construction in Climbing and Descending Tall Isogeny Volcanos (ePrint 2024/924, Sections 3-4) apply unchanged to the degree-73 conductor step over characteristic two? The input is the surface curve y^2 + xy = x^3 + 1 over F_(2^37), registry ID icv1-f2m37-tm534059-32aad96b, trace -534059.",
         x: 45, top: y, width: 522) - 17
y = text("The exact parity gate", x: 45, top: y, width: 522,
         size: 13, bold: true, color: accent) - 8
y = text("The paper's Section 3 requires the Kani degree M and both input isogeny degrees N and m to be coprime to the characteristic. Section 4 sets M = N + m. For odd N, these two remaining requirements cannot hold together in characteristic two:",
         x: 45, top: y, width: 522) - 14

let boxesTop = y
box(CGRect(x: 45, y: boxesTop - 102, width: 249, height: 102),
    title: "Choose odd M", body: "Then m = M - N is even. The auxiliary isogeny violates the stated coprimality condition.")
box(CGRect(x: 318, y: boxesTop - 102, width: 249, height: 102),
    title: "Choose odd m", body: "Then M = N + m is even. Full rank-two etale E[M] is unavailable on an ordinary curve in characteristic two.")
y = boxesTop - 120
y = text("Both branches leave the cited proof's hypotheses. This does not rule out a different characteristic-two construction using a non-etale kernel or another auxiliary map.",
         x: 45, top: y, width: 522) - 16
y = text("Exact C37 calculation", x: 45, top: y, width: 522,
         size: 13, bold: true, color: accent) - 8
y = text("A native Rust check computes pi = tau^37 with tau^2 + tau + 2 = 0. For M = 81 = 9 x 9, m = M - 73 = 8 = 2^2 + 2^2. The candidate's auxiliary degree is even. The modular Frobenius values are:",
         x: 45, top: y, width: 522) - 7
y = text("mod 9: pi = -tau, order 24.     mod 81: pi = 27 - tau, order 216.",
         x: 53, top: y, width: 506, size: 10.5, bold: true) - 9
y = text("Full K0[9] is rational over F_(2^(37 x 24)). This arithmetic neither constructs a degree-73 isogeny nor transfers a point.",
         x: 45, top: y, width: 522)
pageFooter(1)
context.endPDFPage()

context.beginPDFPage(nil)
y = 745
y = text("Source audit and remaining obligation", x: 45, top: y,
         width: 522, size: 20, bold: true, color: accent) - 9
line(y)
y -= 19
y = text("Correction to the earlier E8 note", x: 45, top: y,
         width: 522, size: 13, bold: true, color: accent) - 8
y = text("The first pass treated full rational E[M] over a single compositum as a required input. Section 4 of ePrint 2024/924 instead works prime by prime with Frobenius eigenbases and avoids materializing that compositum in its efficient implementation. Torsion-field costs require an audit of any actual characteristic-two adaptation; the blanket full-M objection was too strong.",
         x: 45, top: y, width: 522) - 18
y = text("Evidence", x: 45, top: y, width: 522,
         size: 13, bold: true, color: accent) - 8
y = text("The protocol and native tool were committed before the retained replay (final source commit 8e7d9a6148df). Three exact-arithmetic tests passed with Rust 1.93.1. The retained JSON receipt is byte-identical to an earlier development calculation; SHA-256: 6118d92f89aef25c0310e83f37405f56f414ed5b394578898ccdbe25f3c83a5e. Source, protocol, JSON and the full report are in research/isogeny_conductor_gap_20261007/.",
         x: 45, top: y, width: 522) - 19
y = text("Remaining work", x: 45, top: y, width: 522,
         size: 13, bold: true, color: accent) - 8
y = text("To advance E8 to the C37 prototype, establish a characteristic-two theorem or algorithm that handles the auxiliary map and kernel, construct an explicit isogeny between registered endpoints, verify its action on the named subgroup, and count complete setup and transfer operations. Only then can the n=131 cost model be evaluated. There is no operation-cost or discrete-logarithm result from this preflight.",
         x: 45, top: y, width: 522) - 19
y = text("Primary source", x: 45, top: y, width: 522,
         size: 13, bold: true, color: accent) - 8
_ = text("Steven D. Galbraith, Climbing and Descending Tall Isogeny Volcanos, extended ePrint 2024/924, Sections 3-4 and Theorem 1. https://eprint.iacr.org/2024/924",
         x: 45, top: y, width: 522)
pageFooter(2)
context.endPDFPage()
context.closePDF()
